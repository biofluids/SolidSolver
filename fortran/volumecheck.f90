module volumecheck
	
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
				end if
			end if
			volume = volume + v(i)
		end do
	end function volume			
		
	function c3d8_det(x,y,z)
		implicit none
		real(8), dimension(3) :: x,y,z
		real(8) :: c3d8_det
		c3d8_det = x(1)*(y(2)*z(3)-y(3)*z(2)) - y(1)*(x(2)*z(3)-x(3)*z(2)) + z(1)*(x(2)*y(3)-x(3)*y(2))
	end function c3d8_det
		
	
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
		volume_c3d8 = c3d8_det(a(:,1),a(:,2),a(:,3)) + c3d8_det(a(:,1),a(:,4),a(:,5)) + c3d8_det(a(:,1),a(:,6),a(:,7))
	end function volume_c3d8
	
end module volumecheck
		