program solidsolver
	
	use read_file
	use integration

	
	implicit none
	
	integer :: reduced = 0, npt
	
	call read_input(10,'input.txt',simu_type, maxit, nsteps, nprint, tol, dt, damp, materialprops, gravity)
	call read_mesh(nsd,nn,nel,nen,coords,connect,bc1,bc2)
	
	write(*,'("simu_type=",i6,5x,"tol=",e6.1,5x,"maxit=",i6,5x)') simu_type,tol,maxit
	write(*,'("nsteps=",i6,5x,"dt=",e6.1,5x,"nprint=",i6,5x,"damp=",f6.2)') nsteps,dt,nprint,damp
	write(*,*) 'materialprops:'
	write(*,'(e10.3)') materialprops(:)
	write(*,*) 'gravity:'
	write(*,'(f6.2)') gravity(:)
	
	write(*,*) nsd, nn, nel, nen, shape(bc1), shape(bc2)
	
	npt = int_number(nsd,nen,reduced)
	
	write(*,*) 'number of integration points', npt
	write(*,*) 'weights of integration', int_weights(nsd, nen, npt)
	write(*,*) 'integration points are: ', int_points(nsd, nen, npt)
	
	
	
	
	
	
end program solidsolver