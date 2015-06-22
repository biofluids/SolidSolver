program test

	implicit none
	
	real(8), dimension(2,3) :: a 
	integer :: i,j
	real(8), dimension(2) :: b
	
	a = reshape([1.,4.,2.,5.,3.,6.],shape(a))
	b = [1.,2.]
	do i=1,3
		a(:,i) = a(:,i) + b
	end do
	write(*,*) a
		
		
		
	
	
end program