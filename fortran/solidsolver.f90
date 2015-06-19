program solidsolver
	
	use read_file
	use integration
	use face
	use shapefunction
	use material
	use externalforce
	use internalforce
	use tangentstiffness
    use mgmres
	
	implicit none
	
	integer :: i,step,nit,row,col,j,k
	real(8), allocatable, dimension(:) :: Fext, F1, F2, Fint, R, w, w1, dw
	real(8), allocatable, dimension(:,:) ::  A
	real(8) :: loadfactor, increment, err1, err2
	
	
	
	
	call read_input(10,'input.txt',simu_type, maxit, firststep, adjust, nsteps, nprint, tol, dt, damp, materialprops, gravity)
	call read_mesh(nsd,ned,nn,nel,nen,coords,connect,bc1,bc2)
	
	allocate(Fext(nn*ned))
	allocate(F1(nn*ned))
	allocate(F2(nn*ned))
	allocate(Fint(nn*ned))
	allocate(R(nn*ned))
	allocate(w(nn*ned))
	allocate(w1(nn*ned))
	allocate(dw(nn*ned))
	allocate(A(nn*ned,nn*ned))
	
	! initialize w
	do i=1,nn*ned
		w(i) = 0.
	end do

	if (simu_type == 0) then
		loadfactor = 0.
		increment = firststep
		step = 1
		
		!call write_results
		call force_traction(F1)
		call force_body(F2)
		Fext = F1 + F2
		
		do while (loadfactor < 1.)
			if (loadfactor+increment>1.) then
				increment = 1. - loadfactor
			end if
			w1 = w
			loadfactor  = loadfactor + increment
			err1 = 1.
			err2 = 1.
			nit = 0
			write(*,'("==============================Step",i5,5x,"Load",f12.4,"==============================")') step,loadfactor
			do while (((err1>tol) .OR. (err2>tol)) .and. (nit<maxit))
				nit = nit + 1
				Fint = force_internal(w)
				A = tangent_internal(w)
				R = Fint - loadfactor*Fext
				! fix the prescribed displacement
				do i=1,size(bc1,2)
					row = ned*int((bc1(1,i)-1.)) + int(bc1(2,i))
					do col=1,ned*nn
						A(row,col) = 0.
						A(col,row) = 0.
					end do
					A(row,col) = 1.
					R(row) = -bc1(3,i) + w(row)
				end do
				! solve
				call solve_mgmres(nn*ned,A,-R,dw)
				w = w + dw
				! check convergence
				err1 = sqrt(dot_product(dw,dw)/dot_product(w,w))
				err2 = sqrt(dot_product(R,R))/(ned*nn)
				write(*,'("Iteration number:",i8,5x,"Err1:",E12.4,5x,"Err2:",E12.4,5x,"Tolerance:",E12.4)') nit,err1,err2,tol
			end do
			step = step + 1
			if (nit == maxit) then
				w = w1
				loadfactor = loadfactor - increment
				increment = increment / 2
			else if (nit > 6) then
				increment = increment*(1-adjust)
			else if (nit < 6) then
				increment = increment*(1+adjust)
			end if
		end do
		
	end if	

	
	
	
	
end program solidsolver