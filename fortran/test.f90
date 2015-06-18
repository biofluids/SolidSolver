program test
	use material
	implicit none
	integer :: nsd,ned
	real(8), dimension(2,2) :: B, stress
	real(8) :: Ja
	real(8), dimension(5) :: materialprops
	real(8), dimension(2,2,2,2) :: C
	
	nsd = 2
	ned = 2
	B = reshape([14,11,15,12],shape(B))
	Ja = 3.
	materialprops = [3.,1.,5.,10.,2.]
	C = materialstiffness(nsd,ned,B,Ja,materialprops)
	stress = Kirchhoffstress(nsd,ned,B,Ja,materialprops)
	
	write(*,*) 'C(:,:,1,1) ', C(:,:,1,1)
	write(*,*) 'C(:,:,2,1) ', C(:,:,2,1)
	write(*,*) 'C(:,:,1,2) ', C(:,:,1,2)
	write(*,*) 'C(:,:,2,2) ', C(:,:,2,2)
	
	write(*,*) 'stress: ', stress

end program