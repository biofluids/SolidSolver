program solidsolver
	
	use read_file
	use integration
	use face
	use shapefunction
	use material
	use externalforce
	use internalforce

	
	implicit none
	
	integer :: reduced = 0, npt,i,row,j
	real(8), allocatable, dimension(:) :: Fext, F1, F2, Fint, dofs
	real(8), allocatable, dimension(:,:) :: temp
	
	call read_input(10,'input.txt',simu_type, maxit, nsteps, nprint, tol, dt, damp, materialprops, gravity)
	call read_mesh(nsd,ned,nn,nel,nen,coords,connect,bc1,bc2)
	
	!write(*,'("simu_type=",i6,5x,"tol=",e6.1,5x,"maxit=",i6,5x)') simu_type,tol,maxit
	!write(*,'("nsteps=",i6,5x,"dt=",e6.1,5x,"nprint=",i6,5x,"damp=",f6.2)') nsteps,dt,nprint,damp
	!write(*,*) 'materialprops:'
	!write(*,'(e10.3)') materialprops(:)
	!write(*,*) 'gravity:'
	!write(*,'(f6.2)') gravity(:)
	
	
	allocate(Fext(nn*ned))
	allocate(F1(nn*ned))
	allocate(F2(nn*ned))
	allocate(Fint(nn*ned))
	allocate(dofs(nn*ned))
	allocate(temp(nsd,nn))
	
	call force_traction(F1)
	call force_body(F2)
	Fext = F1 + F2
	!open(20,file='externalforce.txt')
	!do i=1,nn*ned
	!	write(20,*) Fext(i)
	!end do
	
	open(40,file='displacement.txt')
	do i=1,3
		do j=1,nn
			row = ned*(i-1) + j
			read(40,*) temp(i,j)
		end do
	end do
	do i=1,nn
		do j=1,3
			row=ned*(i-1) + j
			dofs(row) = temp(j,i)
		end do
	end do
			
	
	!Fint = force_internal(dofs)
	!open(30,file='internalforce.txt')
	!do i=1,nn*ned
	!	write(30,*) Fint(i)
	!end do
	
	
	
	
end program solidsolver