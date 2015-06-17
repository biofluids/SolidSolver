module ttt
	implicit none
	
contains
	subroutine xyz(a)
		implicit none
		integer, intent(inout) :: a
		integer :: b
		b = 100
		a = a + b
	end subroutine xyz
	
	subroutine zzz(a)
		implicit none
		integer, intent(inout) :: a
		integer :: b
		write(*,*) b
		a = a + b
	end subroutine zzz
end module ttt