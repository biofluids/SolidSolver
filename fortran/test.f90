program test

	implicit none
	
	character(50):: filepath
	character(50) :: filename
	filepath = '/Users/Jie/Documents/SolidResults/'
	filename = trim(filepath)//'solid.geo'
	write(*,*) filename
	

end program