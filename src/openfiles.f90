!***********************************************************************
!
subroutine openfiles
!
!***********************************************************************
!
  use title_mod
  implicit none
!
!***********************************************************************
!

  ! Simulation log file, monitor file for residuals
  open(unit=6,file=monitor_file)
  rewind 6

end subroutine
