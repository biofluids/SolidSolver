module output
    implicit none
contains
    subroutine write_results(filepath,dofs)
        use read_file, only: mode, step, nn, nsd, nel, isbinary
        character(80) :: filepath
        real(8), dimension(nn*nsd), intent(in) :: dofs
        if (step == 0) then
            call write_case(filepath)
        end if
        call write_geometry(filepath,dofs)
        call write_displacement(filepath,dofs)
        call write_stress(filepath)
    end subroutine write_results

    subroutine write_case(filepath)
        use read_file, only: nsteps, nprint, dt
        character(80) :: filepath, filename
        real(8), dimension(nsteps/nprint+1) :: time
        integer :: i
        
        filename=trim(filepath)//'solid.case'
        open(unit=10,file=trim(filename))
        write(10,'("FORMAT",/)') 
        write(10,'("type:",12x,"ensight gold",/)') 
        write(10,'("GEOMETRY",/)')
        write(10,'("model:",12x,"solid.geo******",12x,"change_coords_only",/)')
        write(10,'("VARIABLE",/)') 
        write(10,'("vector per node:",12x,"displacement",12x,"solid.dis******")')
        write(10,'("tensor symm per node:",12x,"stress",12x,"solid.sig******",/)')
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

    subroutine write_geometry(filepath,dofs)
        use read_file, only: step, nn, nsd, nen, nel, connect, coords, nprint, isbinary
        real(8), dimension(nn*nsd), intent(in) :: dofs
        character(80) :: filepath, filename, buffer
        integer :: i,j
        integer, dimension(nn) :: nodeid
        integer, dimension(nel) :: eleid
        real(8), dimension(3,nn) :: coords1
        character(6) :: temp
        
        write(temp,'(i6.6)') step/nprint
        filename = trim(filepath)//'solid.geo'//trim(temp)
        
        do i=1,nsd
            do j=1,nn
                coords1(i,j) = coords(i,j) + dofs((j-1)*nsd+i)
            end do
        end do
        if (nsd==2) then
            do i=1,nn
                coords1(3,i) = 0.
            end do
        end if
        
        if (isbinary == 1) then
            open(unit=11,file=trim(filename),form='unformatted')
            buffer = 'Fortran Binary'
            write(11) buffer
            buffer = 'This is a geometry file'
            write(11) buffer
            buffer = 'Useless discription line'
            write(11) buffer
            buffer = 'node id given'
            write(11) buffer
            buffer = 'element id given'
            write(11) buffer
            buffer = 'part'
            write(11) buffer
            i = 1
            write(11) i
            buffer = 'solid'
            write(11) buffer
            buffer = 'coordinates'
            write(11) buffer
            write(11) nn
            do i = 1, nn
                nodeid(i) = i
            end do
            write(11) nodeid
            write(11) sngl(coords1(1,:))
            write(11) sngl(coords1(2,:))
            write(11) sngl(coords1(3,:))
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
            write(11) buffer
            write(11) nel
            do i = 1, nel
                eleid(i) = i
            end do
            write(11) eleid
            write(11) connect
            close(11)
        else
            open(unit=10,file=trim(filename),form='FORMATTED')
            write(10,'(A)') 'Ensight Model Geometry File'
            write(10,'(A)') 'Ensight Model Geometry File'
            write(10,'(A)') 'node id given'
            write(10,'(A)') 'element id given'
            write(10,'(A)') 'part'
            write(10,'(i10)') 1
            write(10,'(A)') 'solid'
            write(10,'(A)') 'coordinates'
            write(10,'(i10)') nn
            do i=1,nn
                write(10,'(i10)') i
            end do
            do i=1,3
                do j=1,nn
                    write(10,'(e12.5)') coords1(i,j)
                end do
            end do
            if (nsd==2) then
                if (nen==3) then
                    write(10,'(A)') 'tria3'
                elseif (nen==4) then
                    write(10,'(A)') 'quad4'
                end if
            elseif (nsd==3) then
                if (nen==4) then
                    write(10,'(A)') 'tetra4'
                elseif (nen==8) then
                    write(10,'(A)') 'hexa8'
                end if
            end if
            write(10,'(i10)') nel
            do i=1,nel
                write(10,'(i10)') i
            end do
            do i=1,nel
                write(10,'(*(i10))') connect(:,i)
            end do
            close(10)
        end if
    end subroutine write_geometry


    subroutine write_displacement(filepath,dofs)
        use read_file, only: step, nsd, nsd, nn, coords, nel, nen, connect, nprint, isbinary
        real(8), dimension(nn*nsd), intent(in) :: dofs
        character(80) :: filepath, filename, buffer
        integer :: i,j,row
        character(6) :: temp
        real(8), dimension(3,nn) :: displacement
        
        do i=1,nsd
            do j=1,nn
                row=(j-1)*nsd+i
                displacement(i,j)=dofs(row)
            end do
        end do
        if (nsd==2) then
            do i=1,nn
                displacement(3,i) = 0.
            end do
        end if
        
        write(temp,'(i6.6)') step/nprint
        filename = trim(filepath)//'solid.dis'//trim(temp)
        
        if (isbinary == 1) then
            open(unit=11,file=trim(filename),form='unformatted')
            buffer = 'This is a vector per node file for displacement'
            write(11) buffer
            buffer = 'part'
            write(11) buffer
            i = 1
            write(11) i
            buffer = 'coordinates'
            write(11) buffer
            do i = 1, 3
                write(11) sngl(displacement(i,:))
            end do
            close(11)
        else
            open(unit=10,file=trim(filename),form='FORMATTED')
            write(10,'(A)') 'This is a vector per node file for displacement'
            write(10,'(A)') 'part'
            write(10,'(i10)') 1
            write(10,'(A)') 'coordinates'
            do i=1,3
                do j=1,nn
                    write(10,'(e12.5)') displacement(i,j)
                end do
            end do
            close(10)
        end if
    end subroutine write_displacement

    subroutine write_stress(filepath)
        use read_file, only: step, nn, nprint, isbinary, nodeStress
        character(80) :: filepath, filename, buffer
        character(6) :: temp
        integer :: i, j

        write(temp,'(i6.6)') step/nprint
        filename = trim(filepath)//'solid.sig'//trim(temp)

        ! write to file
        if (isbinary == 1) then
            open(unit=11,file=trim(filename),form='unformatted')
            buffer = 'This is a symm tensor per node file for stress'
            write(11) buffer
            buffer = 'part'
            write(11) buffer
            i = 1
            write(11) i
            buffer = 'coordinates'
            write(11) buffer
            do i = 1, 6
                write(11) sngl(nodeStress(i,:))
            end do
            close(11)
        else
            open(unit=10,file=trim(filename),form='FORMATTED')
            write(10,'(A)') 'This is a symm tensor per node file for stress'
            write(10,'(A)') 'part'
            write(10,'(i10)') 1
            write(10,'(A)') 'coordinates'
            do i=1,6
                do j=1,nn
                    write(10,'(e12.5)') nodeStress(i,j)
                end do
            end do
            close(10)
        end if
    end subroutine write_stress
end module output
