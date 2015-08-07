module read_file
	
	implicit none
	
	integer :: simu_type, maxit, nsteps, nprint, step, isbinary
	real(8) :: firststep, adjust, tol, dt, damp
	real(8) :: materialprops(5), gravity(3)
	integer :: nsd, nen, nn, nel, no_bc1, no_bc2, ned
	integer, allocatable :: connect(:,:), bc1(:,:)
	real(8), allocatable :: coords(:,:), bc2(:,:)
	integer, allocatable :: share(:)
	
	save
	
contains
	subroutine read_input(unitnum,filename,simu_type, maxit, firststep, adjust, nsteps, nprint, tol, dt, damp, materialprops, &
						  gravity, isbinary)
		implicit none
		
		integer, intent(in) :: unitnum
		character(len=*), intent(in) :: filename
	
		character(100) :: text
		character(1) :: flag=':'
		integer :: i, j, k, l, ios
		real(8) :: temp(18)
	
		integer, intent(out) :: simu_type, maxit, nsteps, nprint, isbinary
		real(8), intent(out) :: firststep, adjust, tol, dt, damp, materialprops(5), gravity(3)
	
		open(unit=unitnum,file=filename)
		i=0
		ios=0
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
		firststep=temp(4)
		adjust=temp(5)
		nsteps=int(temp(6))
		dt=temp(7)
		nprint=int(temp(8))
		damp=temp(9)
		materialprops(:)=temp(10:14)
		gravity(:)=temp(15:17)
	    isbinary = temp(18)
		close(10)
	end subroutine read_input
	
	subroutine read_mesh(nsd,ned,nn,nel,nen,coords,connect,bc1,bc2,share)
		implicit none
		
		integer, intent(out) :: nsd, ned, nen, nn, nel
		integer, allocatable, intent(out) :: connect(:,:), bc1(:,:)
		real(8), allocatable, intent(out) :: coords(:,:), bc2(:,:)
		integer, allocatable :: share(:)
		integer :: no_bc1, no_bc2, i,j
	
		open(10,file='coords.txt')
		read(10,*) nsd, nn
		allocate(coords(3,nn)) 
		do i=1,nn
			read(10,*) coords(:,i)
		end do
		close(10)
	
		open(10,file='connect.txt')
		read(10,*) nel, nen
		allocate(connect(nen,nel))
		do i=1,nel
			read(10,*) connect(:,i)
		end do
		close(10)
	
		open(10,file='bc.txt')
		read(10,*) no_bc1
		if (no_bc1 /= 0) then
			allocate(bc1(2,no_bc1))
			do i=1,no_bc1
				read(10,*) bc1(:,i)
			end do
		else
			allocate(bc1(2,1))
			bc1 = -1
		end if	
		close(10)
	
		open(10,file='load.txt')
		read(10,*) no_bc2
		allocate(bc2(5,no_bc2))
		do i=1,no_bc2
			read(10,*) bc2(:,i)
		end do
		
		ned = nsd
		
		allocate(share(nn))
		do i=1,nn
			share(i) = 0
		end do
		do i=1,nel
			do j=1,nen
				share(connect(j,i)) = share(connect(j,i)) + 1
			end do
		end do
		
	    close(10)
	
	end subroutine read_mesh
	
end module read_file
	
	