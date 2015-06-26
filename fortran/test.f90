program test
	use volumecheck
	implicit none
	
	real(8), dimension(3,8) :: coords
	real(8) :: f
	
	coords=reshape([0.,1.,1.,0.,0.,1.,1.,0.,0.,0.,0.,0.,1.,1.,1.,1.,1.,1.,0.,0.,1.,1.,0.,0.],shape(coords))
	f = volume_c3d8(coords)
	write(*,*) f
	
		
		
	
	
end program