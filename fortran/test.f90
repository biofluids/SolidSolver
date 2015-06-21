program test

	implicit none
	
	integer, dimension(2,3) :: a 
	integer :: i,j
	integer, dimension(3) :: b
	
	a = reshape([1,4,2,5,3,6],shape(a))
	do i=1,2
		
			write(*,'(*(i10))') a(i,:)
		
	end do
	
end program