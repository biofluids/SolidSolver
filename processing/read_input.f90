program read_input
	implicit none
	
	character(100) :: text
	character(1) :: flag=':'
	integer :: i, j, k, l, ios=0
	real(8) :: temp(16)
	
	integer simu_type, maxit, nsteps, nprint
	real(8) :: tol, relax, dt, damp
	real(8) :: materialprops(5), gravity(3)
	
	open(10,file='input.txt')
	i=0
	do while(ios==0)
		read(10,'(a)',IOSTAT=ios) text
		j=index(text,flag)
		l=0
		if (j /= 0) then
			i=i+1
			do k=j+1,len_trim(text)
				if (text(k:k) /= ' ') then
					l=l+1
				endif
			end do
			read(text(j+1:j+1+l),*) temp(i)
		end if	
	end do
	
	simu_type=int(temp(1))
	tol=temp(2)
	maxit=int(temp(3))
	relax=temp(4)
	nsteps=int(temp(5))
	dt=temp(6)
	nprint=int(temp(7))
	damp=temp(8)
	materialprops(:)=temp(9:13)
	gravity(:)=temp(14:16)
	
	close(10)
	
	write(*,'("simu_type=",i6,5x,"tol=",e6.1,5x,"maxit=",i6,5x,"relax=",f3.1)') simu_type,tol,maxit,relax
	write(*,'("nsteps=",i6,5x,"dt=",e6.1,5x,"nprint=",i6,5x,"damp=",f6.2)') nsteps,dt,nprint,damp
	
	write(*,*) 'materialprops:'
	write(*,'(e10.3)') materialprops(:)
	write(*,*) 'gravity:'
	write(*,'(f6.2)') gravity(:)
	
end program read_input