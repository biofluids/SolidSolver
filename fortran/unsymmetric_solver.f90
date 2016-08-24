! DEPRECATED
module unsymmetric_solver
	implicit none
contains
    subroutine ma41ds(ndofs, Ain, R, dw)
        implicit none
        integer, intent(in) :: ndofs
        real(8), dimension(ndofs, ndofs), intent(in) :: Ain
    	real(8), dimension(ndofs), intent(in) :: R
    	real(8), dimension(ndofs), intent(out) :: dw
    	! ma41 parameters
    	real(8), dimension(10) :: cntl
    	integer, dimension(20) :: icntl, info
    	integer, dimension(50) :: keep
    	integer :: job, n, ne, maxis, maxs
    	integer, allocatable :: irn(:), jcn(:), is(:)
    	real(8), dimension(20) :: rinfo
    	real(8), allocatable :: aspk(:), rhs(:), colsca(:), rowsca(:), s(:)
    	! trivial parameters
    	integer :: i, j, k
    
    	external ma41id, ma41ad

    	! Initialization of control parameters
    	call MA41ID(cntl, icntl, keep)
    	! Request error analysis and iterative refinement
    	icntl(10) = 10
    	icntl(11) = 1
    
    	! Count the nonzeros
    	n = ndofs
    	ne = 0
    	k = 0
    	do j = 1, n
        	do i = 1, n
            	if (Ain(i,j) /= 0.) then
                	ne = ne + 1
            	end if
        	end do
    	end do
    	allocate(irn(ne))
    	allocate(jcn(ne))
    	allocate(aspk(ne))
    	allocate(rhs(n))
    	rhs = R
    	do j = 1, n
        	do i = 1, n
            	if (Ain(i,j) /= 0.) then
                	k = k + 1
                	irn(k) = i
                	jcn(k) = j
                	aspk(k) = Ain(i,j)
            	end if
        	end do
    	end do
    	allocate(colsca(n))
    	allocate(rowsca(n))
   
    	! Analyse sparsity pattern using minimum degree ordering
    	! Scale (column scaling) and factorize matrix
    	! Solve the equations, and perform error analysis
    	job = 6
		maxis = 10*max(info(7), 2*ne + 11*n + 1)
		maxs = 2*maxis
		allocate(is(maxis))
		allocate(s(maxs))
    	call ma41ad(job, n, ne, irn, jcn, aspk, rhs, colsca, rowsca, keep, is, maxis, s, maxs, cntl, icntl, info, rinfo)
		dw = rhs

    	deallocate(aspk)
    	deallocate(rhs)
    	deallocate(colsca)
    	deallocate(rowsca)
    	deallocate(s)
    	deallocate(irn)
    	deallocate(jcn)
    	deallocate(is)

    end subroutine ma41ds
	
end module unsymmetric_solver
