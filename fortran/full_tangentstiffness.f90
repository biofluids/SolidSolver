module full_tangentstiffness

	implicit none
	
contains
	function tangent_internal(dofs)
		use read_file, only: nsd, ned, nn, coords, nel, nen, connect, materialprops
		use shapefunction
		use integration
		use full_material
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
		real(8), allocatable, dimension(:,:) :: xilist
		real(8), allocatable, dimension(:) :: weights
		integer :: ele,a,i,npt,j,row,intpt,l,d,k,col
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
		tangent_internal = 0.
		! allocate
		npt = int_number(nsd,nen,0)
		allocate(xilist(nsd,npt))
		allocate(weights(npt))
		xilist = int_points(nsd,nen,npt)
		weights = int_weights(nsd,nen,npt)
				
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
			kint = 0.
			! loop over integration points
			do intpt=1,npt
				xi = xilist(:,intpt)
				dNdxi = sfder(nen,nsd,xi)
				! set up the jacobian matrix
				dxdxi = matmul(elecoord,dNdxi)
				! inverse matrix and determinant
				dxidx = dxdxi
				call DGETRF(n1,n1,dxidx,n1,ipiv,info)
				call DGETRI(n1,dxidx,n1,ipiv,work,n1,info)
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
				call DGETRI(n1,Finv,n1,ipiv,work,n1,info)
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
											            (eye(i,k)*stress(j,l)+C(i,j,k,l))*dNdy(a,j)*dNdy(d,l)*weights(intpt)*det;			
									end do
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
		
		deallocate(xilist)
		deallocate(weights)
		
	end function tangent_internal
end module full_tangentstiffness
		