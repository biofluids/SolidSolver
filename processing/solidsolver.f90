program solidsolver
	
	use read_file
	
	implicit none
	
	call read_input(10,'input.txt',simu_type, maxit, nsteps, nprint, tol, dt, damp, materialprops, gravity)
	call read_mesh(nsd,nn,nel,nen,coords,connect,bc1,bc2)
	
	write(*,'("simu_type=",i6,5x,"tol=",e6.1,5x,"maxit=",i6,5x)') simu_type,tol,maxit
	write(*,'("nsteps=",i6,5x,"dt=",e6.1,5x,"nprint=",i6,5x,"damp=",f6.2)') nsteps,dt,nprint,damp
	write(*,*) 'materialprops:'
	write(*,'(e10.3)') materialprops(:)
	write(*,*) 'gravity:'
	write(*,'(f6.2)') gravity(:)
	
	write(*,*) nsd, nn, nel, nen, shape(bc1), shape(bc2)
	
	
	
	
end program solidsolver