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
  use k_epsilon_std
  use temperature, only: t
  use statistics

  implicit none
!
!***********************************************************************
!

!--------------------------------------------------------------
!    [ Writing restart file]
!--------------------------------------------------------------
  open(unit=3,file=restart_file,form='unformatted')
  rewind 3

  write(3) itime,time
  if(const_mflux) write(3) gradpcmf
  write(3) flmass
  write(3) u
  write(3) v
  write(3) w
  write(3) p
  write(3) te
  write(3) ed
  write(3) t
  write(3) vis
  !      write(3) vart
  !      write(3) edd
  !      write(3) ret
  !      write(3) den
  !      write(3) utt
  !      write(3) vtt
  !      write(3) wtt
  write(3) uu
  write(3) vv
  write(3) ww
  write(3) uv
  write(3) uw
  write(3) vw
  write(3) uo
  write(3) vo
  write(3) wo
  !      write(3) to
  write(3) teo
  write(3) edo
  !      write(3) varto
  !      write(3) con
  !      write(3) cono
  !      write(3) alph

  rewind 3
  close (3)

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

  write(6,*)'=*=*= Simulation restart files have been written! =*=*='

end subroutine
