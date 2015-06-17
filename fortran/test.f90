program test
	integer, dimension(3,2) :: a
	
	a=reshape([1,1,1,2,2,2],shape(a))
	write(*,*) size(a,2)
end program