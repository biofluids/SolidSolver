else
	allocate(M(nn*ned,nn*ned))
	allocate(Kint(nn*ned,nn*ned))
	allocate(eye(nn*ned,nn*ned))
	allocate(un(nn*ned))
	allocate(vn(nn*ned))
	allocate(un1(nn*ned))
	allocate(vn1(nn*ned))
	allocate(an(nn*ned))
	allocate(an1(nn*ned))
	allocate(F(nn*ned))
	gamma = 0.5
	beta = 0.25
	! initialize
	do i=1,nn*ned
		un(i) = 0.
		un1(i) = 0.
		vn(i) = 0.
		vn1(i) = 0.
		an(i) = 0.
		an1(i) = 0.
	end do
	call write_results(filepath,w)
	! mass matrix and external force are constant
	call mass_matrix(M)
	call force_traction(F1)
	call force_body(F2)
	Fext = F1 + F2
	! initial Fint
	Fint = force_internal(w)
	F = Fext - Fint
	! in order to get an right, modify
	do i=1,size(bc1,2)
		row = ned*int((bc1(1,i)-1.)) + int(bc1(2,i))
		do col=1,nn*ned
			M(row,col) = 0.
			M(col,row) = 0.
		end do
		M(row,row) = 1.
		F(row) = 0.
	end do
	! initial acceleration
	call solve_mgmres(nn*ned,M,F,an)
	do step = 1,nsteps
		write(*,'("==============================Step",i5,5x,"Time",f12.4,"==============================")') step,step*dt
		err1 = 1.
		err2 = 1.
		nit = 0
		w = un
		do while( ((err1>tol) .or. (err2>tol)) .and. (nit<maxit)  )
			nit = nit + 1
			an1 = (w - (un+dt*vn+0.5*dt**2*(1-2*beta)*an))/(beta*dt**2)
			vn1 = vn + (1-gamma)*dt*an + gamma*dt*an1
			Fint = force_internal(w)
			Kint = tangent_internal(w)
			F = Fext - Fint - damp*vn1
			R = matmul(M,an1) - F
			! modify to get A and R right
			do i=1,size(bc1,2)
				row = ned*int((bc1(1,i)-1.)) + int(bc1(2,i))
				do col=1,nn*ned
					Kint(row,col) = 0.
					Kint(col,row) = 0.
				end do
				Kint(row,row) = 1.-1./(beta*dt**2)
				R(row) = w(row) - bc1(3,i)
			end do
			! Jacobian
			do i=1,nn*ned
				do j=1,nn*ned
					if (i == j) then
						eye(i,j)=1.
					else
						eye(i,j)=0.
					end if
				end do
			end do
			A = (M+damp*gamma*dt*eye)/(beta*dt**2) + Kint
			call solve_mgmres(nn*ned,A,-R,dw)
			! check convergence
			w = w + dw
			err1 = sqrt(dot_product(dw,dw)/dot_product(w,w))
			err2 = sqrt(dot_product(R,R))/(ned*nn)
		end do
		write(*,'("Iteration number:",i8,5x,"Err1:",E12.4,5x,"Err2:",E12.4,5x,"Tolerance:",E12.4)') nit,err1,err2,tol
		vn = vn1
		un = w
		an = an1
		
		if (MOD(step,nprint)==0) then
			call write_results(filepath,w)
		end if
	end do
end if	
