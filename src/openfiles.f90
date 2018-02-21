!***********************************************************************
!
subroutine openfiles
!
!***********************************************************************
!
  use title_mod
  implicit none
  
  ! character(8)  :: date
  ! character(10) :: time
!
!***********************************************************************
!

  ! Simulation log file

  open(unit=6,file=monitor_file)
  rewind 6


  ! Open folder with data for postprocessing in Paraview
  
  ! call date_and_time(DATE=date, TIME=time)
  ! write(datetime, '(a)') date(1:4)//"-"//date(5:6)//"-"//date(7:8)//"_"//time(1:2)//":"//time(3:4)//":"//time(5:6)

  call execute_command_line("mkdir VTK")


end subroutine
