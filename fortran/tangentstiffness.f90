module tangentstiffness
	
	
	implicit none
	
contains
	function tangent_internal(dofs)
		use read_file, only: nsd, ned, nn, coords, nel, nen, connect, materialprops
		use shapefunction
		use integration
		use material
		implicit none
		real(8), dimension(nn*ned), intent(in) :: dofs
		real(8), dimension(nn*ned,nn*ned) :: tangent_internal
		real(8), dimension(ned,nsd,ned,nsd) :: C
		real(8), dimension(nsd,nen) :: elecoord
		real(8), dimension(ned,nen) :: eledof
		real(8), dimension(nen,nsd) :: dNdx, dNdy
		real(8), dimension(ned,nsd) :: stress
		real(8), dimension(nen*ned,nen*ned) :: kint
		real(8), dimension(nsd) :: xi
		real(8), dimension(nen,nsd) :: dNdxi 
		real(8), dimension(nsd,nsd) :: dxdxi, dxidx, F, Finv, B, eye
		real(8), allocatable, dimension(:,:) :: xilist, xilist1
		real(8), allocatable, dimension(:) :: weights, weights1
		integer :: ele,a,i,npt,j,row,intpt,npt1,l,d,k,col
		real(8) :: det, Ja
		real(8), dimension(nsd) :: work ! for lapack inverse
		integer, dimension(nsd) :: ipiv ! for lapack inverse
		integer :: info, n1 ! for lapack inverse
		
		external DGETRF
		external DGETRI
		n1 = nsd
		! square matrix
		do i=1,nsd
			do j=1,nsd
				if(i==j) then
					eye(i,j) = 1.
				else
					eye(i,j) = 0.
				end if
			end do
		end do
		! initialize 
		do i=1,ned*nn
			do j=1,ned*nn
				tangent_internal(i,j) = 0.
			end do
		end do
		! allocate
		if (.NOT. allocated(xilist)) then
			npt = int_number(nsd,nen,0)
			allocate(xilist(nsd,npt))
			allocate(weights(npt))
		end if	
		! loop over elements
		do ele=1,nel
			! extract coords of nodes, and dofs for the element
			do a=1,nen
				elecoord(:,a) = coords(:,connect(a,ele))
				do i=1,ned
					eledof(i,a) = dofs(ned*(connect(a,ele)-1)+i)
				end do
			end do
			! initialize
			do i=1,nen*ned
				do j=1,nen*ned
					kint(i,j) = 0.
				end do
			end do
			! fully integration
			! set up integration points and weights
			npt = int_number(nsd,nen,0)
			xilist = int_points(nsd,nen,npt)
			weights = int_weights(nsd,nen,npt)
			! loop over integration points
			do intpt=1,npt
				xi = xilist(:,intpt)
				dNdxi = sfder(nen,nsd,xi)
				! set up the jacobian matrix
				dxdxi = matmul(elecoord,dNdxi)
				! inverse matrix and determinant
				dxidx = dxdxi
				call DGETRF(n1,n1,dxidx,n1,ipiv,info)
				if (info /= 0) then
					stop 'Matrix is numerically singular!'
				end if
				call DGETRI(n1,dxidx,n1,ipiv,work,n1,info)
				if (info /= 0) then
					stop 'Matrix inversion failed!'
				end if	
				!
				if (nsd == 2) then
					det = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
				else if (nsd == 3) then
					det = dxdxi(1,1)*dxdxi(2,2)*dxdxi(3,3) - dxdxi(1,1)*dxdxi(3,2)*dxdxi(2,3) &
						  - dxdxi(1,2)*dxdxi(2,1)*dxdxi(3,3) + dxdxi(1,2)*dxdxi(2,3)*dxdxi(3,1) &
						  + dxdxi(1,3)*dxdxi(2,1)*dxdxi(3,2) - dxdxi(1,3)*dxdxi(2,2)*dxdxi(3,1)
				end if
				! compute dNdx
				dNdx = matmul(dNdxi,dxidx)
				! deformation gradient, F_ij = delta_ij + dU_i/dx_j
				F = eye + matmul(eledof,dNdx)
				! left Cauchy-Green tensor, B = FF^T and Ja = det(F)
				B = matmul(F,transpose(F))
				if (nsd == 2) then
					Ja = F(1,1)*F(2,2) - F(1,2)*F(2,1)
				else if (nsd == 3) then
					Ja = F(1,1)*F(2,2)*F(3,3) - F(1,1)*F(3,2)*F(2,3) &
						  - F(1,2)*F(2,1)*F(3,3) + F(1,2)*F(2,3)*F(3,1) &
						  + F(1,3)*F(2,1)*F(3,2) - F(1,3)*F(2,2)*F(3,1)
				end if
				! compute dNdy, in which y is the coord. after deformation
				! inverse of F
				Finv = F
				call DGETRF(n1,n1,Finv,n1,ipiv,info)
				if (info /= 0) then
					stop 'Matrix is numerically singular!'
				end if
				call DGETRI(n1,Finv,n1,ipiv,work,n1,info)
				if (info /= 0) then
					stop 'Matrix inversion failed!'
				end if	
				!
				dNdy = matmul(dNdx,Finv)
				! compute the Kirchhoff stress
				stress = Kirchhoffstress(nsd,ned,B,Ja,materialprops)
				! compute the material stiffness C_ijkl
				C = materialstiffness(nsd,ned,B,Ja,materialprops)
				! compute the element internal force
				do a=1,nen
					do i=1,ned
						do d=1,nen
							do k=1,ned
								row = (a-1)*ned + i
								col = (d-1)*ned + k
								do j=1,nsd
									do l=1,nsd
										kint(row,col) = kint(row,col) + &
											            C(i,j,k,l)*dNdy(a,j)*dNdy(d,l)*weights(intpt)*det;
										kint(row,col) = kint(row,col) - &
										                C(j,j,k,l)/nsd*dNdy(a,i)*dNdy(d,l)*weights(intpt)*det;				
									end do
									kint(row,col) = kint(row,col) - &
									                stress(i,j)*dNdy(a,k)*dNdy(d,j)*det*weights(intpt);
									kint(row,col) = kint(row,col) + &
									                stress(j,j)/nsd*dNdy(a,k)*dNdy(d,i)*det*weights(intpt);
								end do
							end do
						end do
					end do
				end do
			end do
			! reduced integration
			! set up the integration points and weights
			! allocate
			if (.NOT. allocated(xilist1)) then
				npt1 = int_number(nsd,nen,1)
				allocate(xilist1(nsd,npt1))
				allocate(weights1(npt1))
			end if
			xilist1 = int_points(nsd,nen,npt1)
			weights1 = int_weights(nsd,nen,npt1)
			do intpt=1,npt1
				xi = xilist1(:,intpt)
				dNdxi = sfder(nen,nsd,xi)
				! jacobian
				dxdxi = matmul(elecoord,dNdxi)
				! inverse matrix and determinant
				dxidx = dxdxi
				call DGETRF(n1,n1,dxidx,n1,ipiv,info)
				if (info /= 0) then
					stop 'Matrix is numerically singular!'
				end if
				call DGETRI(n1,dxidx,n1,ipiv,work,n1,info)
				if (info /= 0) then
					stop 'Matrix inversion failed!'
				end if	
				!
				if (nsd == 2) then
					det = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
				else if (nsd == 3) then
					det = dxdxi(1,1)*dxdxi(2,2)*dxdxi(3,3) - dxdxi(1,1)*dxdxi(3,2)*dxdxi(2,3) &
						  - dxdxi(1,2)*dxdxi(2,1)*dxdxi(3,3) + dxdxi(1,2)*dxdxi(2,3)*dxdxi(3,1) &
						  + dxdxi(1,3)*dxdxi(2,1)*dxdxi(3,2) - dxdxi(1,3)*dxdxi(2,2)*dxdxi(3,1)
				end if
				! compute dNdx
				dNdx = matmul(dNdxi,dxidx)
				! deformation gradient, F_ij = delta_ij + dU_i/dx_j
				F = eye + matmul(eledof,dNdx)
				! left Cauchy-Green tensor, B = FF^T and Ja = det(F)
				B = matmul(F,transpose(F))
				if (nsd == 2) then
					Ja = F(1,1)*F(2,2) - F(1,2)*F(2,1)
				else if (nsd == 3) then
					Ja = F(1,1)*F(2,2)*F(3,3) - F(1,1)*F(3,2)*F(2,3) &
						  - F(1,2)*F(2,1)*F(3,3) + F(1,2)*F(2,3)*F(3,1) &
						  + F(1,3)*F(2,1)*F(3,2) - F(1,3)*F(2,2)*F(3,1)
				end if
				! compute dNdy, in which y is the coord. after deformation
				! inverse of F
				Finv = F
				call DGETRF(n1,n1,Finv,n1,ipiv,info)
				if (info /= 0) then
					stop 'Matrix is numerically singular!'
				end if
				call DGETRI(n1,Finv,n1,ipiv,work,n1,info)
				if (info /= 0) then
					stop 'Matrix inversion failed!'
				end if			
				dNdy = matmul(dNdx,Finv)
				! compute the Kirchhoff stress
				stress = Kirchhoffstress(nsd,ned,B,Ja,materialprops)
				! compute the material stiffness C_ijkl
				C = materialstiffness(nsd,ned,B,Ja,materialprops)
				! compute the element internal force
				do a=1,nen
					do i=1,ned
						do d=1,nen
							do k=1,ned
								row = (a-1)*ned + i
								col = (d-1)*ned + k
								do j=1,nsd
									do l=1,nsd
										kint(row,col) = kint(row,col) + &
										                C(j,j,k,l)/nsd*dNdy(a,i)*dNdy(d,l)*weights1(intpt)*det;				
									end do
									kint(row,col) = kint(row,col) - &
									                stress(j,j)/nsd*dNdy(a,k)*dNdy(d,i)*det*weights1(intpt);
								end do
							end do
						end do
					end do
				end do
			end do
			! scatter the element kint into the global Kint
			!
			do a=1,nen
				do i=1,ned
					do d=1,nen
						do k=1,ned
							row = ned*(connect(a,ele)-1) + i;
							col = ned*(connect(d,ele)-1) + k;
		                    tangent_internal(row,col) = tangent_internal(row,col) + kint(ned*(a-1)+i,ned*(d-1)+k);
						end do
					end do
				end do
			end do
		end do
	end function tangent_internal
end module tangentstiffness
		