program test
	!use volumecheck
	implicit none
	
	real(8), dimension(3,8) :: coords
	real(8) :: f
	
	coords=reshape([0.,0.,0.,1.,0.,0.,1.,1.,0.,0.,1.,0.,0.,0.,1.,1.,0.,2.,1.,1.,2.,0.,1.,2.],shape(coords))
	f = volume_c3d8(coords)
	write(*,*) f
	f = volume_c3d8_new(coords)
	write(*,*) f
	
		
		
	
	
contains


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
	real(8), dimension(3,9) :: A
	A(:,1) = x(:,7) - x(:,1)
	A(:,2) = x(:,2) - x(:,1)
	A(:,3) = x(:,3) - x(:,6)
	A(:,4) = A(:,1)
	A(:,5) = x(:,5) - x(:,1)
	A(:,6) = x(:,6) - x(:,8)
	A(:,7) = A(:,1)
	A(:,8) = x(:,4) - x(:,1)
	A(:,9) = x(:,8) - x(:,3)
	volume_c3d8=(det3(A(:,1),A(:,2),A(:,3))+det3(A(:,4),A(:,5),A(:,6))+det3(A(:,7),A(:,8),A(:,9)))/6.
end function volume_c3d8

function volume_c3d8_new(x)
	implicit none
	real(8), dimension(3,8), intent(in) :: x
	real(8) :: volume_c3d8_new
	real(8), dimension(3,6) :: A
	A(:,1) = x(:,7) - x(:,1)
	A(:,2) = x(:,2) + x(:,8) - x(:,5) - x(:,6)
	A(:,3) = x(:,8) - x(:,2) + x(:,5) - x(:,6)
	A(:,4) = A(:,1)
	A(:,5) = x(:,2) + x(:,8) - x(:,3) - x(:,4)
	A(:,6) = x(:,2) - x(:,8) + x(:,3) - x(:,4)
	volume_c3d8_new=(det3(A(:,1),A(:,2),A(:,3))+det3(A(:,4),A(:,5),A(:,6)))/12.
end function volume_c3d8_new

end program