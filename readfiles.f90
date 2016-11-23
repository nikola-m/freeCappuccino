!***********************************************************************
!
subroutine readfiles
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use variables
  use title_mod
  use statistics
  use k_epsilon_std, only: te,ed,teo,edo 
  use temperature, only: t

  implicit none  
!
!***********************************************************************
!

  open(unit=3,file=restart_file,form='unformatted')
  rewind 3

  read(3) itime,time
  if(const_mflux) read(3) gradpcmf
  read(3) flmass(:)
  read(3) u(:)
  read(3) v(:)
  read(3) w(:)
  read(3) p(:)
  read(3) te(:)
  read(3) ed(:)
  read(3) t(:)
  read(3) vis(:)
  !read(3) vart(:)
  !read(3) edd(:)
  !read(3) ret(:)
  !read(3) den(:)
  !read(3) utt(:)
  !read(3) vtt(:)
  !read(3) wtt(:)
  read(3) uu(:)
  read(3) vv(:)
  read(3) ww(:)
  read(3) uv(:)
  read(3) uw(:)
  read(3) vw(:)
  read(3) uo(:)
  read(3) vo(:)
  read(3) wo(:)
  !read(3) to(:)
  read(3) teo(:)
  read(3) edo(:)
  !read(3) varto(:)
  !read(3) con(:)
  !read(3) cono(:)
  !read(3) alph(:)

  rewind 3
  
  close (3)

  if (ltransient) then
  !------------------------------------------------
  !     [read statistics after first collection: ]
  !------------------------------------------------
  !      open(unit=85,file=trim(out_folder_path)//'/statistics1')   ! <- n_sample is here, statistics restart file 1
  !      open(unit=86,file=trim(out_folder_path)//'/statistics2')   ! <- u_aver, v_aver,... are here, statistics restart file 2
  !      rewind 85
  !      rewind 86

  !      read(85,*) n_sample
  !      read(86,*) u_aver,v_aver,w_aver, &
  !                 uu_aver,vv_aver,ww_aver, &
  !                 uv_aver,uw_aver,vw_aver, &
  !                 te_aver
  !      close (85)
  !      close (86)
  endif

end subroutine
