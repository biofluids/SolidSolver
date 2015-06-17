module external
	implicit none
	
contains
	subroutine externalforce(Fext)
		! returns the external force
		use read_file, only: nsd,ned,nn,nel,nen,coords,connect,bc1,bc2, gravity, materialprops
		use face
		use shapefunction
		use integration
		
		implicit none
		
		real(8), dimension(nn*ned), intent(out) :: Fext
		
		integer :: j,i, no_bc2, ele, faceid, nfacenodes, a, npt, intpt, row
		integer, allocatable, dimension(:) :: nodelist,nodelist1
		real(8), allocatable, dimension(:,:) :: coord, xilist, dNdxi, xilist1
		real(8), allocatable, dimension(:) :: f, weights, N, weights1
		real(8), dimension(nsd) :: traction, xi
		real(8), dimension(ned*nen) :: fb
		real(8), dimension(nen) :: N1
		real(8), dimension(nsd,nsd) :: dxdxi
		real(8), dimension(nsd,nen) :: elecoord
		real(8), dimension(nen,nsd) :: dNdxi1
		real(8) :: det, rho
		
		no_bc2 = size(bc2,2)
		! allocate
		if (.NOT. allocated(nodelist)) then
			nfacenodes = face_nodes_no(nsd,nen)
			allocate(nodelist(nfacenodes))
			allocate(coord(nsd,nfacenodes))
			allocate(f(ned*nfacenodes))
			npt = int_number(nsd-1, nfacenodes, 0)
			allocate(xilist(nsd,npt))
			allocate(weights(npt))
			allocate(dNdxi(nfacenodes,nsd-1))
			allocate(N(nfacenodes))
		end if
		
		!
		! force due to traction
		! loop over faces with prescribed traction
		do j=1,no_bc2
			! extract the coords and traction of that face
			ele = int(bc2(1,j))
			faceid = int(bc2(2,j))
			nfacenodes = face_nodes_no(nsd,nen)
			nodelist = face_nodes(nsd,nen,nfacenodes,faceid)
			do a=1,nfacenodes
				coord(:,a) = coords(:,connect(nodelist(a),ele))
			end do
			do i=1,nsd
				traction(i) = bc2(i+2,j)
			end do
			! compute the force on the face
			do i=1,ned*nfacenodes
				f(i) = 0.
			end do
			! set up integration points and weights
			npt = int_number(nsd-1, nfacenodes, 0)
			xilist = int_points(nsd-1,nfacenodes,npt)
			weights = int_weights(nsd-1,nfacenodes,npt)
			
			! loop over the integration points
			do intpt=1,npt
				xi(:) = xilist(:,intpt)
				dNdxi(:,:) = sfder(nfacenodes,nsd-1,xi)
				N(:) = sf(nfacenodes,nsd-1,xi)
				! set up the jacobian matrix
				dxdxi(:,:) = matmul(coord,dNdxi)
				if (nsd == 2) then
					det = sqrt((dxdxi(1,1))**2 + (dxdxi(2,1)**2))
				else if (nsd == 3) then
					det = sqrt( ((dxdxi(2,1)*dxdxi(3,2))-(dxdxi(2,2)*dxdxi(3,1)))**2 &
                		  + ((dxdxi(1,1)*dxdxi(3,2))-(dxdxi(1,2)*dxdxi(3,1)))**2 &
               		      + ((dxdxi(1,1)*dxdxi(2,2))-(dxdxi(1,2)*dxdxi(2,1)))**2)
				end if
				! compute the force
				do a=1,nfacenodes
					do i=1,nsd
						row = ned*(a-1)+i
						f(row) = f(row) + traction(i)*N(a)*weights(intpt)*det
					end do
				end do
			end do
			! scatter the force into the global external force
			do a=1,nfacenodes
				do i=1,ned
					row = (connect(nodelist(a),ele)-1)*ned + i
					Fext(row) = Fext(row) + f((a-1)*ned+i)
				end do
			end do
		end do
		
		
		if (.NOT. allocated(xilist1)) then
			npt = int_number(nsd, nen, 0)
			allocate(xilist1(nsd,npt))
			allocate(weights1(npt))
		end if
		!
		! force due to gravity
		! loop over elements
		do ele=1,nel
			! extract coords of the nodes, and dof for the element
			do a=1,nen
				elecoord(:,a) = coords(:,connect(a,ele))
			end do
			! compute the force on this element
			do i=1,ned*nen
				fb(i) = 0.
			end do
			rho = materialprops(5)
			! set up integration points and weights
			npt = int_number(nsd, nen, 0)
			xilist1 = int_points(nsd,nen,npt)
			weights1 = int_weights(nsd,nen,npt)
			! loop over integration points
			do intpt=1,npt
				xi(:) = xilist1(:,intpt)
				dNdxi1 = sfder(nen,nsd,xi)
				N1 = sf(nen,nsd,xi)
				! set up the jacobian matrix
				dxdxi(:,:) = matmul(elecoord,dNdxi)
				det = dxdxi(1,1)*dxdxi(2,2) - dxdxi(1,2)*dxdxi(2,1)
				! compute the force
				do a=1,nen
					do i=1,ned
						row=ned*(a-1) + i
						fb(row) = fb(row) + gravity(i)*rho*N1(a)*weights1(intpt)*det
					end do
				end do
				! scatter
				do a=1,nen
					do i=1,ned
						row = ned*(connect(a,ele)-1)+i;
						Fext(row) = Fext(row) + fb(ned*(a-1)+i);
					end do
				end do
			end do
		end do
			
	end subroutine externalforce
end module external
					  
				
			