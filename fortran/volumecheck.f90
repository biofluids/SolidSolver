module volumecheck
	
	! This module was originally designed to compute the volume of the whole part before and after deformation.
	! But I could not find a nice formulation to calculate the volume of an arbitrary hexahedron given the
	! coordinates of its nodes.
	! Therefore the volume change rate of individual elements are checked in the output subroutine instead,
	! assuming all the elements are of the same size, average volume change rate is computed.
	! This module is not in use until the volume formulation is found.  
    ! Jie Cheng, 08/13/15
	
	implicit none
	
contains
	
		
	function volume(w)
		use read_file, only: nsd, nen, ned, coords, connect, nn, nel
		implicit none
		real(8), dimension(nn*nsd), intent(in) :: w
		real(8) :: volume
		real(8), dimension(nel) :: v
		real(8), dimension(nsd,nn) :: coords1
		real(8), dimension(nsd,nen) :: elecoord
		integer :: i,j
		
		volume = 0.
		
		do i=1,nsd
			do j=1,nn
				coords1(i,j) = coords(i,j) + w((j-1)*nsd+i)
			end do
		end do
		
		
		do i=1,nel
			do j=1,nen
				elecoord(:,j) = coords1(:,connect(j,i))
			end do
			if (nsd==3) then
				if (nen==8) then
					v(i) = volume_c3d8(elecoord)
				else if (nen==4) then
					v(i) = volume_c3d4(elecoord)
				end if
			else if (nsd==2) then
				if (nen==4) then
					v(i) = volume_c2d4(elecoord)
				else if (nen==3) then
					v(i) = volume_c2d3(elecoord)
				end if	
			end if
			volume = volume + v(i)
		end do
	end function volume			
		
	function det3(x,y,z)
		implicit none
		real(8), dimension(3) :: x,y,z
		real(8) :: det3
		det3 = x(1)*(y(2)*z(3)-y(3)*z(2)) - y(1)*(x(2)*z(3)-x(3)*z(2)) + z(1)*(x(2)*y(3)-x(3)*y(2))
	end function det3
		
	
	function volume_c3d8(x)
		implicit none
		real(8), dimension(3,8), intent(in) :: x
		real(8) :: volume_c3d8
		real(8), dimension(3,7) :: a
		a(:,1) = x(:,7) - x(:,1)
		a(:,2) = x(:,2) - x(:,1)
		a(:,3) = x(:,3) - x(:,6)
		a(:,4) = x(:,5) - x(:,1)
		a(:,5) = x(:,6) - x(:,8)
		a(:,6) = x(:,4) - x(:,1)
		a(:,7) = x(:,8) - x(:,3)
		volume_c3d8 = (det3(a(:,1),a(:,2),a(:,3)) + det3(a(:,1),a(:,4),a(:,5)) + det3(a(:,1),a(:,6),a(:,7)))/6.
	end function volume_c3d8
	
	function volume_c3d4(x)
		implicit none
		real(8), dimension(3,4), intent(in) :: x
		real(8) :: volume_c3d4
		real(8), dimension(3,3) :: a
		a(1,1) = x(1,1) - x(1,4)
		a(2,1) = x(2,1) - x(2,4)
		a(3,1) = x(3,1) - x(3,4)
		a(1,2) = x(1,2) - x(1,4)
		a(2,2) = x(2,2) - x(2,4)
		a(3,2) = x(3,2) - x(3,4)
		a(1,3) = x(1,3) - x(1,4)
		a(2,3) = x(2,3) - x(2,4)
		a(3,3) = x(3,3) - x(3,4)
		volume_c3d4 = det3(a(:,1),a(:,2),a(:,3))/6.
	end function volume_c3d4
	
	function volume_c2d4(x)
		implicit none
		real(8), dimension(2,4), intent(in) :: x
		real(8) :: volume_c2d4
		real(8), dimension(2,2) :: a
		a(:,1)=x(:,3)-x(:,1)
		a(:,2)=x(:,4)-x(:,2)
		volume_c2d4=0.5*abs(a(1,1)*a(2,2)-a(1,2)*a(2,1))
	end function
	
	function volume_c2d3(x)
		implicit none
		real(8), dimension(2,3), intent(in) :: x
		real(8) :: volume_c2d3
		volume_c2d3 = 0.5*abs(x(1,1)*(x(2,2)-x(2,3)) + x(1,2)*(x(2,3)-x(2,1)) + x(1,3)*(x(2,1)-x(2,2)))
	end function
		
	
end module volumecheck
		