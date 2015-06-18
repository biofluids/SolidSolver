module internalforce
	
	
	implicit none
	
contains
	function force_internal(dofs)
		use read_file, only: nsd, ned, nn, coords, nel, nen, connect, materialprops
		use shapefunction
		use integration
		use material
		implicit none
		real(8), dimension(nn*ned), intent(in) :: dofs
		real(8), dimension(nn*ned) :: Fint
		real(8), dimension(nsd,nen) :: elecoord
		real(8), dimension(ned,nen) :: eledof
		real(8), dimension(nen,nsd) :: dNdx, dNdy
		real(8), dimension(ned,nsd) :: stress
		real(8), dimension(nen*ned) :: fele
		real(8), dimension(nsd) :: xi
		real(8), dimension(nen,nsd) :: dNdxi 
		real(8), dimension(nsd,nsd) :: dxdxi, dxidx, F, Finv, B, eye
		real(8), allocatable, dimension(:,:) :: xilist
		real(8), allocatable, dimension(:) :: weights
		integer :: ele,a,i,npt,j,row,intpt
		real(8) :: det, Ja
		
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
		
		! initialize Fint
		do i=1,nn*ned
			Fint(i) = 0.
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
			! compute the internal force
			! initialize
			do i=1,ned*nen
				fele(i) = 0.
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
				do i=1,nsd
					do j=1,nsd
						dxidx(i,j) = dxdxi(j,i)
					end do 
				end do
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
					Ja = B(1,1)*B(2,2) - B(1,2)*B(2,1)
				else if (nsd == 3) then
					Ja = B(1,1)*B(2,2)*B(3,3) - B(1,1)*B(3,2)*B(2,3) &
						  - B(1,2)*B(2,1)*B(3,3) + B(1,2)*B(2,3)*B(3,1) &
						  + B(1,3)*B(2,1)*B(3,2) - B(1,3)*B(2,2)*B(3,1)
				end if
				! compute dNdy, in which y is the coord. after deformation
				! inverse of F
				do i=1,nsd
					do j=1,nsd
						Finv(i,j) = F(j,i)
					end do
				end do
				dNdy = matmul(dNdx,Finv)
				! compute the Kirchhoff stress
				stress = Kirchhoffstress(ned,nsd,B,Ja,materialprops)
				! compute the element internal force
				do a=1,nen
					do i=1,nsd
						row=(a-1)*ned+i
						do j=1,nsd
							fele(row) = fele(row) + stress(i,j)*dNdy(a,j)*weights(intpt)*det
							fint(row) = fint(row) - stress(j,j)/nsd*dNdy(a,i)*weights(intpt)*det;
						end do
					end do
				end do
			end do
			!
			!
			!
			!
		end do
	end function force_internal
end module internalforce
				
				
				
				
				
				
				
				
				
				
				
				
				
	