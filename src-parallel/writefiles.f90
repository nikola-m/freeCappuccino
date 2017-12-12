!***********************************************************************
!
subroutine writefiles
!
!***********************************************************************
!
  use types
  use parameters
  use title_mod
  use geometry
  use variables
  use statistics
  use sparse_matrix
  use output

  implicit none
!
!***********************************************************************
!
  integer :: i
  integer :: output_unit
  character( len = 5) :: nproc_char

  ! nproc_char <- myid zapisan levo u vidu stringa.
  call i4_to_s_left ( myid, nproc_char )

  ! Write in a char variable the current timestep number and create a folder with this name
  write(timechar,'(i6)') itime
  call execute_command_line("mkdir processor"//trim(nproc_char)//"/"//adjustl(timechar))
  

  ! Open an output file in this folder for velocity
  call get_unit ( output_unit )
  open ( unit = output_unit, file = trim(out_folder_path)//"/processor"//trim(nproc_char)//'/'//trim(adjustl(timechar))//'/U')

  rewind output_unit

  ! Write output to file
  do i=1,numCells
    write(output_unit,'(a,es11.4,1x,es11.4,1x,es11.4,a)') '(',u(i),v(i),w(i),')'
  enddo

  close(output_unit)


  ! Write output files in Paraview .vtu format
  !+-----------------------------------------------------------------------------+
  call get_unit( output_unit )

  open(unit = output_unit, file = trim(out_folder_path)//"/processor"//trim(nproc_char)// &
&                                 '/'//trim(adjustl(timechar))//'/velocity_field.vtu')

  call vtu_write_vector_field( output_unit, 'velocity_field', u, v, w )

  close (  unit = output_unit )



  ! Open an output file in this folder for pressure
  call get_unit ( output_unit )
  open ( unit = output_unit, file = trim(out_folder_path)//"/processor"//trim(nproc_char)//'/'//trim(adjustl(timechar))//'/p')

  rewind output_unit

  ! Write output to file
  do i=1,numCells
    write(output_unit ,'(es11.4)') p(i)
  enddo

  close(output_unit)


  ! Write output files in Paraview .vtu format
  !+-----------------------------------------------------------------------------+
  call get_unit( output_unit )

  open(unit = output_unit, file = trim(out_folder_path)//"/processor"//trim(nproc_char)// &
&                                 '/'//trim(adjustl(timechar))//'/pressure_field.vtu')

  call vtu_write( output_unit, 'pressure_field', p )

  close (  unit = output_unit )


  if(solveTKE) then

    ! Open an output file in this folder for Turbulence kinetic energy
    call get_unit ( output_unit )
    open ( unit = output_unit, file = trim(out_folder_path)//"/processor"//trim(nproc_char)//'/'//trim(adjustl(timechar))//'/k')

    rewind output_unit

    ! Write output to file
    do i=1,numCells
      write(output_unit ,'(es11.4)') te(i)
    enddo

    close(output_unit)


    ! Write output files in Paraview .vtu format
    !+-----------------------------------------------------------------------------+
      call get_unit( output_unit )

      open(unit = output_unit, file = trim(out_folder_path)//"/processor"//trim(nproc_char)// &
    &                                 '/'//trim(adjustl(timechar))//'/tke_field.vtu')

      call vtu_write( output_unit, 'tke_field', te )

      close (  unit = output_unit )


  endif



  if( solveEpsilon ) then

    ! Open an output file in this folder for Epsilon
    call get_unit ( output_unit )
    open ( unit = output_unit, file = trim(out_folder_path)//"/processor"//trim(nproc_char)// &
  &                                   '/'//trim(adjustl(timechar))//'/epsilon')

    rewind output_unit

    ! Write output to file
    do i=1,numCells
      write(output_unit ,'(es11.4)') ed(i)
    enddo

    close(output_unit)


    ! Write output files in Paraview .vtu format
    !+-----------------------------------------------------------------------------+
      call get_unit( output_unit )

      open(unit = output_unit, file = trim(out_folder_path)//"/processor"//trim(nproc_char)// &
    &                                 '/'//trim(adjustl(timechar))//'/epsilon_field.vtu')

      call vtu_write( output_unit, 'epsilon_field', ed )

      close (  unit = output_unit )


  endif

  

  if( solveOmega ) then

    ! Open an output file in this folder for Omega
    call get_unit ( output_unit )
    open ( unit = output_unit, file = trim(out_folder_path)//"/processor"//trim(nproc_char)//'/'//trim(adjustl(timechar))//'/omega')

    rewind output_unit

    ! Write output to file
    do i=1,numCells
      write(output_unit ,'(es11.4)') ed(i)
    enddo

    close(output_unit)


    ! Write output files in Paraview .vtu format
    !+-----------------------------------------------------------------------------+
      call get_unit( output_unit )

      open(unit = output_unit, file = trim(out_folder_path)//"/processor"//trim(nproc_char)// &
    &                                 '/'//trim(adjustl(timechar))//'/omega_field.vtu')

      call vtu_write( output_unit, 'omega_field', ed )

      close (  unit = output_unit )


  endif      

end subroutine
