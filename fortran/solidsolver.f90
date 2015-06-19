program solidsolver
	
	use read_file
	use integration
	use face
	use shapefunction
	use material
	use externalforce
	use internalforce
	use tangentstiffness

	
	implicit none
	
	integer :: i,j
	real(8), allocatable, dimension(:) :: Fext, F1, F2, Fint, dofs
	real(8), allocatable, dimension(:,:) ::  Kint
	
	call read_input(10,'input.txt',simu_type, maxit, nsteps, nprint, tol, dt, damp, materialprops, gravity)
	call read_mesh(nsd,ned,nn,nel,nen,coords,connect,bc1,bc2)
	
	allocate(Fext(nn*ned))
	allocate(F1(nn*ned))
	allocate(F2(nn*ned))
	allocate(Fint(nn*ned))
	allocate(dofs(nn*ned))
	allocate(Kint(nn*ned,nn*ned))
	
	call force_traction(F1)
	call force_body(F2)
	Fext = F1 + F2

	
	
	
	
end program solidsolver