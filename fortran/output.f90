module output
	implicit none
	
contains
	subroutine write_results(dofs)
		use read_file, only: simu_type, step, nn, ned
		implicit none
		real(8), dimension(nn*ned), intent(in) :: dofs
		if (step == 0) then
			call write_case()
		end if
		call write_geometry(dofs)
		
	end subroutine write_results
	subroutine write_case()
		use read_file, only: nsteps, nprint, dt, ned, nn
		implicit none
		
		
		real(8), dimension(nsteps/nprint+1) :: time
		integer :: i
		
		
		open(10,file='/Users/jiecheng/Documents/SolidResults/solid.case')
		write(10,'("FORMAT",/)') 
		write(10,'("type:",12x,"ensight gold",/)') 
		write(10,'("GEOMETRY",/)')
		write(10,'("model:",12x,"solid.geo******",12x,"change_coords_only",/)')
		!write(10,'("VARIABLE",/)') 
		!write(10,'("vector per node:",12x,"displacement",12x,"solid.dis******")')
		!write(10,'("tensor symm per node:",12x,"stress",12x,"solid.sig******",/)')
		write(10,'("TIME",/)')
		write(10,'("time set:",12x,i10)') 1
		write(10,'("number of steps:",12x,i10)') nsteps/nprint + 1
		write(10,'("filename start number:",12x,i10)') 0
		write(10,'("filename increment:",12x,i10)') 1
		write(10,'("time values:")')
		do i=1,(nsteps/nprint+1)
			time(i) = (i-1)*nprint*dt
			write(10,'(f12.3)')  time(i)
		end do
		close(10)
	end subroutine write_case
	
	subroutine write_geometry(dofs)
		use read_file, only: step, nn, nsd,ned, nen, nel, connect, coords
		implicit none
		
		real(8), dimension(nn*ned), intent(in) :: dofs
		integer :: i,j
		real(8), dimension(3,nn) :: coords1
		character, dimension(80) :: buffer
		integer, dimension(nn) :: node_id
		integer, dimension(nel) :: element_id
		
		character(54) :: filename, temp
		
		write(temp,'(i6.6)') step
		filename = '/Users/jiecheng/Documents/SolidResults/solid.geo'//trim(temp)
		
		
		do i=1,nn
			node_id(i) = i
		end do
		do i=1,nel
			element_id = i
		end do
		do i=1,nsd
			do j=1,nn
				coords1(i,j) = coords(i,j) + dofs((j-1)*ned+i)
			end do
		end do
		if (nsd==2) then
			do i=1,nn
				coords1(3,i) = 0.
			end do
		end if
		i = 1
		open(10,file=filename,form='UNFORMATTED')
		buffer = 'Fortran Binary'
		write(10) buffer
		buffer = 'Ensight Model Geometry File'
		write(10) buffer
		buffer = 'Ensight Model Geometry File'
		write(10) buffer
		buffer = 'node id given'
		write(10) buffer
		buffer = 'element id given'
		write(10) buffer
		buffer = 'part'
		write(10) buffer
		write(10,'(i10)') i
		buffer = 'solid'
		write(10) buffer
		buffer = 'coordinates'
		write(10) buffer
		write(10,'(i10)') nn
		write(10,'(i10)') node_id
		do i=1,3
			write(10,'(E12.5)') sngl(coords1(i,:))
		end do
		if (nsd==2) then
		    if (nen==3) then
				buffer = 'tria3'
		    elseif (nen==4) then
				buffer = 'quad4'
		    end if
		elseif (nsd==3) then
		    if (nen==4) then
				buffer = 'tetra4'
		    elseif (nen==8) then
				buffer = 'hexa8'
		    end if
		end if
		write(10) buffer
		write(10,'(i10)') nel
		write(10,'(i10)') element_id
		do i=1,nel
			write(10,'(i10)') connect(:,i)
		end do
	end subroutine write_geometry
	
end module output
			
			
		
		
		
		
		
		