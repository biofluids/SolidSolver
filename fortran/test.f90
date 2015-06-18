program test
	use material
	implicit none
	integer, dimension(3,3) :: eye
	integer, dimension(3,3) :: a, b, c
	integer :: i,j
	do i=1,3
		do j=1,3
			if (i==j) then
				eye(i,j) = 1
			else 
				eye(i,j) = 0
			end if
		end do
	end do
	
	a = reshape([1,2,3,4,5,6,7,8,9],shape(a))
	!b = a+eye
	
	c = matmul(a,transpose(a))
	
	write(*,*) c
	

end program