!***********************************************************************
!
subroutine write_restart_files
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use variables
  use title_mod
  use statistics
  use utils

  implicit none

  integer :: restart_unit
  character( len = 5) :: nproc_char 

!
!***********************************************************************
!

  ! NOTE: nproc_char <- this (=myid + 1) written as left aligned string.
  call i4_to_s_left ( myid, nproc_char )

  call get_unit ( restart_unit )

  open ( unit = restart_unit,file=adjustl(trim(restart_file))//'-'//trim(nproc_char),form='unformatted')
  rewind restart_unit

  write(restart_unit) itime,time
  if(const_mflux) write(restart_unit) gradpcmf
  write(restart_unit) flmass
  write(restart_unit) u
  write(restart_unit) v
  write(restart_unit) w
  write(restart_unit) p
  write(restart_unit) te
  write(restart_unit) ed
  write(restart_unit) t
  write(restart_unit) vis
  !      write(restart_unit) vart
  !      write(restart_unit) edd
  !      write(restart_unit) ret
  !      write(restart_unit) den
  !      write(restart_unit) utt
  !      write(restart_unit) vtt
  !      write(restart_unit) wtt
  write(restart_unit) uu
  write(restart_unit) vv
  write(restart_unit) ww
  write(restart_unit) uv
  write(restart_unit) uw
  write(restart_unit) vw
  write(restart_unit) uo
  write(restart_unit) vo
  write(restart_unit) wo
  !      write(restart_unit) to
  write(restart_unit) teo
  write(restart_unit) edo
  !      write(restart_unit) varto
  !      write(restart_unit) con
  !      write(restart_unit) cono
  !      write(restart_unit) alph

  rewind restart_unit
  close (restart_unit)

  if (ltransient) then
!--------------------------------------------------------------
!    [ Writing of the statistics restart file]
!--------------------------------------------------------------
    open(unit=85,file=trim(out_folder_path)//'/statistics1')   ! <- n_sample is here, statistics restart file 1
    open(unit=86,file=trim(out_folder_path)//'/statistics2')   ! <- u_aver, v_aver,... are here, statistics restart file 2
    rewind 85
    rewind 86

    write(85,*) n_sample
    write(86,*) u_aver,v_aver,w_aver, &
                uu_aver,vv_aver,ww_aver, &
                uv_aver,uw_aver,vw_aver,te_aver, &
                te_aver
    close (85)
    close (86)

  endif

  if( myid.eq.0 ) write(6,*)'=*=*= Simulation restart files have been written. =*=*='

end subroutine
