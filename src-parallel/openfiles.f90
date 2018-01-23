!***********************************************************************
!
subroutine openfiles
!
!***********************************************************************
!
  use title_mod
  use parameters, only: nproc
  use utils, only: i4_to_s_left
  implicit none

  character(8)  :: date
  character(10) :: time
  character( len = 5) :: nproc_char
  integer :: id

!
!***********************************************************************
!

  ! Simulation log file, monitor file for residuals
  open(unit=6,file=monitor_file)
  rewind 6



  ! Open folder with data for postprocessing in Paraview
  call date_and_time(DATE=date, TIME=time)

  write(datetime, '(a)') date(1:4)//"-"//date(5:6)//"-"//date(7:8)//"_"//time(1:2)//":"//time(3:4)//":"//time(5:6)

  call execute_command_line("mkdir VTK_"//datetime)

  ! Create folders for process data
  do id=0,nproc-1

    ! nproc_char <- myid zapisan levo u vidu stringa.
    call i4_to_s_left ( id, nproc_char )

    call execute_command_line("mkdir VTK_"//datetime//'/processor'//trim(nproc_char) )
 
  enddo

end subroutine
