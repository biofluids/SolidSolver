program test
	implicit none
	! these should be passed in
	integer :: nsd=3, nen=8
	!
	character(100) :: text
	character(len=:), allocatable :: temp
	integer :: begin_node=0, end_node, begin_element, end_element, nset=0, ios=0, nn, nel, increment
	real(8), allocatable :: coords(:,:)
	integer, allocatable :: connect(:,:)
	
	integer :: line,i,j,k,l
	
	type set
	character(20) name
	integer :: nodes(2000)
	integer :: elements(2000)
	integer :: face2ele(2,3000)
    end type
	
	type(set) boundary(10)
	
	! find out the beginning and endding line numbers of nodes and elements
	open(10,file='hexa_test.inp')
	do
		read(10,'(a)') text
		begin_node = begin_node+1
		if (index(text,'*') == 0) then
			 exit 
		end if
	end do
	end_node = begin_node
	do
		read(10,'(a)') text
		end_node=end_node+1
		if (index(text,'*') /= 0) then
			 exit 
		end if
	end do
	end_node = end_node - 1
	nn = end_node - begin_node + 1
	
	begin_element = end_node+2
	end_element = begin_element
	do
		read(10,'(a)') text
		end_element = end_element+1
		if (index(text,'*') /= 0) then
			 exit 
		end if
	end do
	end_element = end_element - 2
	nel = end_element - begin_element + 1
	
	! find out the number of sets
	if (index(text,'*Nset') /= 0) then
		nset = nset + 1	
	end if
	do while (ios==0)
		read(10,'(a)',IOSTAT=ios) text
		if (index(text,'*Nset') /= 0) then
			nset = nset + 1	
		end if	
	end do
	
	write(*,'("nn =",i8,/,"nel =",i7,/,"nset =",i6)') nn, nel, nset
	close(10)
	
	! 2nd time to read this file
	open(10,file='hexa_test.inp')
	allocate(coords(nn,3))
	allocate(connect(nel,nen))
	do line=1,begin_node-1
		read(10,'(a)') text
	end do
	do i=1,nn
		read(10,*) line,coords(i,:)
	end do
	read(10,'(a)') text
	do i=1,nel
		read(10,*) line,connect(i,:)
		!write(*,*) connect(i,:)
	end do
	! read sets
	do line=1,1
		read(10,'(a)') text
		i=index(text,'=')+1
		text=text(i:len(text))
		j=index(text,',')-1
		boundary(line)%name(:)= text(1:j)
		if (index(text,'generate') /= 0) then
			read(10,*) begin_node,end_node,increment
			boundary(line)%nodes=0.0_8
			do i=1,(end_node-begin_node)/increment+1
				boundary(line)%nodes(i)=begin_node+increment*(i-1)
			end do
		end if
		
		write(*,'(a)') boundary(line)%name
		write(*,*) boundary(line)%nodes(1:9)
		
	end do
	
end program test