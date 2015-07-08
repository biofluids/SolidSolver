program solidsolver
	use read_file
	
	implicit none
	character(80) :: filepath
	integer :: ct, ct_rate, ct_max, ct1  
	
	call timestamp()
	
	filepath = '/Users/Jie/Documents/SolidResults/'
	call system_clock(ct,ct_rate,ct_max)
	
	call read_input(10,'input.txt',simu_type, maxit, firststep, adjust, nsteps, nprint, tol, dt, damp, materialprops, gravity)
	call read_mesh(nsd,ned,nn,nel,nen,coords,connect,bc1,bc2,share)
	
	if (simu_type == 0) then 
		call statics(filepath)
	end if
	
	call system_clock(ct1)
	call timestamp()
	write(*,'("time elapsed:",f12.2,3x,"seconds")') , dble(ct1-ct)/dble(ct_rate)
	
end program solidsolver
	
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
	use volumecheck
	
	implicit none
	
	integer :: i,nit,row,col,j,k, gamma, beta
	real(8), allocatable, dimension(:) :: Fext, F1, F2, Fint, R, F, w, w1, dw
	real(8), allocatable, dimension(:,:) ::  A
	real(8) :: loadfactor, increment, err1, err2
	real(8) :: v, v1
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
	v = volume(w)
	
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
			! fix the prescribed displacement
			do i=1,size(bc1,2)
				row = ned*((bc1(1,i)-1.)) + (bc1(2,i))
				do col=1,ned*nn
					A(row,col) = 0.
					A(col,row) = 0.
				end do
				A(row,row) = 1.
				R(row) = w(row)
			end do
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
	v1 = volume(w)
	write(*,*) '================================================================================================' 
	write(*,'("Total volume change:",e12.4)') v1/v - 1.
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
	
	integer :: i,nit,row,col,j,k, gamma, beta
	real(8), allocatable, dimension(:) :: Fext, F1, F2, Fint, R, F, w, w1, dw
	real(8), allocatable, dimension(:,:) ::  A, m
	real(8) :: loadfactor, increment, err1, err2
	real(8) :: v, v1
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
	allocate(m(nn*ned,nn*ned))
	
	call mass_matrix(m)
	
	! initialize w
	w = 0.
	v = volume(w)
	
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
			! fix the prescribed displacement
			do i=1,size(bc1,2)
				row = ned*((bc1(1,i)-1.)) + (bc1(2,i))
				do col=1,ned*nn
					A(row,col) = 0.
					A(col,row) = 0.
				end do
				A(row,row) = 1.
				R(row) = w(row)
			end do
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
	v1 = volume(w)
	write(*,*) '================================================================================================' 
	write(*,'("Total volume change:",e12.4)') v1/v - 1.
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