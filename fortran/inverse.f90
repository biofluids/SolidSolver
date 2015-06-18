program inverse
	
	implicit none
	
	real(8), dimension(3,3) :: A
	real(8), dimension(3,3) :: Ainv
	
	real(8), dimension(3) :: work 
	integer, dimension(3) :: ipiv
	integer :: n, info
	
	external DGETRF
	external DGETRI
	
	A = reshape([9.,7.,4.,2.,7.,2.,3.,9.,6.],shape(A))
	
	Ainv = A
	n = size(A,1)
	
	call DGETRF(n,n,Ainv,n,ipiv,info)
	
	if (info /= 0) then
		stop 'matrix is numerically singular!'
	end if
	
	call DGETRI(n,Ainv,n,ipiv,work,n,info)
	
	if (info /= 0) then
		stop 'Matrix inversion failed!'
	end if
	
	!write(*,*) Ainv
	
	A = matmul(A,Ainv)
	write(*,*) A
end program inverse
	