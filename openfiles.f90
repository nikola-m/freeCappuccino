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


  open(unit=66,file=monitor_file)
  rewind 66

  ! open(unit=80,file=trim(out_folder_path)//'/U_END')
  ! rewind 80

  ! open(unit=81,file=trim(out_folder_path)//'/PLOTFILE')
  ! rewind 81

  ! open(unit=90,file='NUSSELT')
  ! rewind 90

end subroutine
