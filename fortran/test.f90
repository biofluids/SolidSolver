program test

	implicit none
	
	integer :: i
	i=0
	do while(i<3)
		write(*,*) 'hello'
		i=i+1
	end do
	write(*,*) i
		
		
	
	
end program