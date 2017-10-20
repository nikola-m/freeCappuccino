program boundary

 use boundary_module

implicit none

integer :: boundary_file
integer :: boundary_of
integer :: process_file
character(len=5) :: ch2

integer :: i,j,np
integer :: narg
character(len=30) :: argument
logical :: lpar = .false.



  ! Check arguments
  narg=command_argument_count()

  i = 1
  j = 1

  arg_loop: do
    if ( i .le. narg) then 
       call get_command_argument(i,argument)
    else
        exit arg_loop
    endif

    select case ( trim(argument) )

          case ('-parallel')
            write(*,'(a)') " Parallel run!"
            lpar = .true.
            i = i+1

          case ('np')
            i = i+1
            call get_command_argument(i,argument)
            read(argument,*) np
            write(*,'(a)') " Number of processes: ", np
            i = i+1

          case ('-d')
            i = i+1
            call get_command_argument(i,key(j))
            i = i + 1
            call get_command_argument(i,value(j))
            i = i + 1
            write(*,'(4a)') " Dictionary: ", trim(key(j)), " --> ", trim(value(j))
            j = j + 1

    end select

  enddo arg_loop

  ! Open files to read
  if(lpar) then
    do i=0,np-1,1
      write(ch2,'(i5)') i
      call execute_command_line("cp processor"//trim(adjustl(ch2))//"/constant/polyMesh/boundary &
      &                             processor"//trim(adjustl(ch2))//"/constant/polyMesh/boundary-OF")
      call get_unit( boundary_of )
      open(unit = boundary_of, file='processor'//trim(adjustl(ch2))//'/constant/polyMesh/boundary-OF')
      rewind boundary_of
      call get_unit( boundary_file )
      open(unit = boundary_file, file='processor'//trim(adjustl(ch2))//'/constant/polyMesh/boundary')
      rewind boundary_file
      write(boundary_file, '(a)') '# bcType nFaces startFace'
      call get_unit( process_file )
      open(unit = process_file, file='processor'//trim(adjustl(ch2))//'/constant/polyMesh/process')
      rewind process_file
      write(process_file, '(a)') '# neighbProcNo nFaces startFace'
      !--------------------------------------------------------
      call process_boundary_file(boundary_of,boundary_file,process_file)
      !--------------------------------------------------------
      close( boundary_file )
    enddo
  else
    call execute_command_line("cp polyMesh/boundary polyMesh/boundary-OF")
    call get_unit( boundary_of )
    open(unit = boundary_of, file='polyMesh/boundary-OF')
    rewind boundary_of
    call get_unit( boundary_file )
    open(unit = boundary_file, file='polyMesh/boundary')
    rewind boundary_file
    write(boundary_file, '(a)') '# bcType nFaces startFace'
    !--------------------------------------------------------
    call process_boundary_file(boundary_of,boundary_file)
    !--------------------------------------------------------
    close( boundary_file )
  endif

end program