program solidsolver
	use read_file
	
	implicit none
	character(80) :: filepath
	integer :: ct, ct_rate, ct_max, ct1  
	
	call timestamp()
	
	filepath = '/Users/jiecheng/Documents/SolidResults/'
	call system_clock(ct,ct_rate,ct_max)
	
	call read_input(10,'input.txt',simu_type, maxit, firststep, adjust, nsteps, nprint, tol, dt, damp, &
					materialprops, gravity, isbinary,penalty)
	call read_mesh(nsd,ned,nn,nel,nen,coords,connect,bc1,bc2,share)

	if (simu_type == 0) then 
		call statics(filepath)
	else if (simu_type == 1) then
		call dynamics(filepath)
	else if (simu_type == 2) then
		call debug(filepath) 
	end if
	
	call system_clock(ct1)
	call timestamp()
	write(*,'("time elapsed:",f12.2,3x,"seconds")') , dble(ct1-ct)/dble(ct_rate)
	
end program solidsolver

subroutine debug(filepath)
	! This is to read displacment file from Abaqus and output the stress	
	use read_file
	use output
!	use volumecheck

	implicit none
	
	integer :: i,nit,row,col,j,k
	real(8), allocatable, dimension(:,:) :: u
	real(8), allocatable, dimension(:) :: w
!	real(8) :: v, v1
	character(80) :: filepath
	
	allocate(u(nsd+1,nn))
	allocate(w(nn*ned))
	
	! initialize w
	w = 0.
!	v = volume(w)
	
	nprint=1
	nsteps=1
	dt=1.
	step = 0
	call write_results(filepath,w)

	! read displacment file
	open(40,file='abaqus.rpt')
	do i = 1, nn
		read(40,*) u(:,i)
	end do
	close(40)
	
	do i = 1, nn
		do j = 1, nsd
			k = nsd*(i-1) + j
			w(k) = u(j+1,i)
		end do
	end do
	
	step = nprint
	call write_results(filepath,w)
!	v1 = volume(w)
!	write(*,'("Total volume change:",e12.4)') v1/v - 1.
	
	deallocate(w)
	deallocate(u)
end subroutine debug
	
subroutine statics(filepath)	
	use read_file
	use integration
	use face
	use shapefunction
	use material
	use externalforce
	use internalforce
	use tangentstiffness
    use mgmres
	use directsolver
	use output
!	use volumecheck
	use full_internalforce
	use full_tangentstiffness
	
	implicit none
	
	integer :: i,nit,row,col,j,k
	real(8), allocatable, dimension(:) :: Fext, F1, F2, Fint, R, w, w1, dw
	real(8), allocatable, dimension(:,:) ::  A
	real(8) :: loadfactor, increment, err1, err2
	real(8), dimension(size(bc1,2)) :: constraint
	real(8), dimension(size(bc1,2),nn*ned) :: der_constraint
!	real(8) :: v, v1
	character(80) :: filepath
	
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
	w = 0.
	constraint = 0.
	der_constraint = 0.
!	v = volume(w)
	
	nprint=1
	nsteps=1
	dt=1.
	loadfactor = 0.
	increment = firststep
	step = 0
	call write_results(filepath,w)
	call force_traction(F1)
	call force_body(F2)
	Fext = F1 + F2
	
	do while (loadfactor < 1.)
		step = step + 1
		if (loadfactor+increment>1.) then
			increment = 1. - loadfactor
		end if
		w1 = w
		loadfactor  = loadfactor + increment
		err1 = 1.
		err2 = 1.
		nit = 0
		write(*,'("==============================Step",i5,5x,"Load",e12.4,"====================================")') step,loadfactor
		do while (((err1>tol) .OR. (err2>tol)) .and. (nit<maxit))
			nit = nit + 1
			Fint = force_internal(w)
			A = tangent_internal(w)
			R = Fint - loadfactor*Fext
			! Apply the bc with penalty method
			if (size(bc1,2) /= 0) then
				do i = 1, size(bc1,2)
					row = ned*(int(bc1(1,i))-1) + int(bc1(2,i))
					constraint(i) = w(row) - bc1(3,i)
					der_constraint(i,row) = 1. 
				end do
				! Explicit multiplication is used to avoid large sparse matrix multiplication
				! Actually
				! A = A + penalty*matmul(transpose(der_constraint),der_constraint)
				! R = R + penalty*matmul(transpose(der_constraint),constraint)
				do i = 1, size(bc1,2)
					row = ned*(int(bc1(1,i))-1) + int(bc1(2,i))
					A(row,row) = A(row,row) + penalty*der_constraint(i,row)*der_constraint(i,row)
					R(row) = R(row) + penalty*der_constraint(i,row)*constraint(i)
				end do
				
			end if
			! solve
			!call solve_mgmres(nn*ned,A,-R,dw)
			call ma57ds(nn*ned,A,-R,dw)
			w = w + dw
			! check convergence
			err1 = sqrt(dot_product(dw,dw)/dot_product(w,w))
			err2 = sqrt(dot_product(R,R))/(ned*nn)
			write(*,'("Iteration number:",i8,5x,"Err1:",E12.4,5x,"Err2:",E12.4,5x,"Tolerance:",E12.4)') nit,err1,err2,tol
		end do
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
	step = nprint
	call write_results(filepath,w)
!	v1 = volume(w)
	write(*,*) '================================================================================================' 
!	write(*,'("Total volume change:",e12.4)') v1/v - 1.
	
	deallocate(Fext)
	deallocate(F1)
	deallocate(F2)
	deallocate(Fint)
	deallocate(R)
	deallocate(w)
	deallocate(w1)
	deallocate(dw)
	deallocate(A)
end subroutine statics

subroutine dynamics(filepath)	
	use read_file
	use integration
	use face
	use shapefunction
	use material
	use externalforce
	use internalforce
	use tangentstiffness
    use mgmres
	use directsolver
	use output
	use mass
	use volumecheck
	
	implicit none
	
	integer :: i,nit,row,col,j,k
	real(8), allocatable, dimension(:) :: Fext, F1, F2, Fint, R, F, w, w1, dw, un, un1, vn, vn1, an, an1
	real(8), allocatable, dimension(:,:) ::  A, M, Kint, eye
	real(8) :: loadfactor, increment, err1, err2
	real(8) :: v, v1, gamma, beta
	character(80) :: filepath
	
	allocate(Fext(nn*ned))
	allocate(F1(nn*ned))
	allocate(F2(nn*ned))
	allocate(Fint(nn*ned))
	allocate(R(nn*ned))
	allocate(w(nn*ned))
	allocate(w1(nn*ned))
	allocate(dw(nn*ned))
	allocate(A(nn*ned,nn*ned))
	allocate(M(nn*ned,nn*ned))
	allocate(Kint(nn*ned,nn*ned))
	allocate(eye(nn*ned,nn*ned))
	allocate(F(nn*ned))
	allocate(un(nn*ned))
	allocate(un1(nn*ned))
	allocate(vn(nn*ned))
	allocate(vn1(nn*ned))
	allocate(an(nn*ned))
	allocate(an1(nn*ned))
	
	call mass_matrix(m)
	
	! initialize 
	w = 0.
	un = 0.
	un1 = 0.
	vn = 0.
	vn1 = 0.
	an = 0.
	an1 = 0.
	gamma = 0.5
	beta = 0.25
	v = volume(w)
	
	call write_results(filepath,w)
	call mass_matrix(M)
	call force_traction(F1)
	call force_body(F2)
	Fext = F1 + F2
	Fint = force_internal(w)
	F = Fext - Fint
	
	! in order to get an right, modify
	do i=1,size(bc1,2)
		row = ned*(int(bc1(1,i))-1) + int(bc1(2,i))
		do col=1,ned*nn
			M(row,col) = 0.
			M(col,row) = 0.
		end do
		M(row,row) = 1.
		F(row) = 0.
	end do
	
	! Initial acceleration
	call ma57ds(nn*ned,M,F,an)
	
	do step = 1,nsteps
		write(*,'("==============================Step",i5,5x,"Time",e12.4,"====================================")') step,step*dt
		err1 = 1.
		err2 = 1.
		nit = 0
		w = un
		do while (((err1>tol) .OR. (err2>tol)) .and. (nit<maxit))
			nit = nit + 1
			an1 = (w - (un+dt*vn+0.5*dt**2*(1-2*beta)*an))/(beta*dt**2)
			vn1 = vn + (1-gamma)*dt*an + gamma*dt*an1
			Fint = force_internal(w)
			Kint = tangent_internal(w)
			F = Fext - Fint - damp*vn1
			R = matmul(M,an1) - F
			! Modify to get A and R right
			do i=1,size(bc1,2)
				row = ned*(int(bc1(1,i))-1) + int(bc1(2,i))
				do col=1,ned*nn
					Kint(row,col) = 0.
					Kint(col,row) = 0.
				end do
				Kint(row,row) = 1. - 1./(beta*dt**2)
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
			call ma57ds(nn*ned,A,-R,dw)
			! check convergence
			w = w + dw
			err1 = sqrt(dot_product(dw,dw)/dot_product(w,w))
			err2 = sqrt(dot_product(R,R))/(ned*nn)
		end do
		write(*,'("Iteration number:",i8,5x,"Err1:",E12.4,5x,"Err2:",E12.4,5x,"Tolerance:",E12.4)') nit,err1,err2,tol
		vn = vn1
		un = w
		an = an1
		v1 = volume(w)
		write(*,'("Total volume change:",e12.4)') v1/v - 1.
		if (MOD(step,nprint)==0) then
			call write_results(filepath,w)
		end if
	end do

	deallocate(Fext)
	deallocate(F1)
	deallocate(F2)
	deallocate(Fint)
	deallocate(R)
	deallocate(w)
	deallocate(w1)
	deallocate(dw)
	deallocate(A)
	deallocate(M)
	deallocate(Kint)
	deallocate(eye)
	deallocate(F)
	deallocate(un)
	deallocate(un1)
	deallocate(vn)
	deallocate(vn1)
	deallocate(an)
	deallocate(an1)
	
end subroutine dynamics

subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
! 
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end subroutine timestamp