!***********************************************************************
!
subroutine read_input_file
!
!***********************************************************************
!
! Open & Read and Process Input File
!
!***********************************************************************
  use types
  use parameters
  use gradients, only: lstsq, lstsq_qr, lstsq_dm, gauss, limiter
  use title_mod

  implicit none

  integer :: i,imon
  character(len=2) :: trpn
  character(len=25) :: convective_scheme
!
!***********************************************************************
!

  OPEN(UNIT=5,FILE=input_file)
  REWIND 5

  READ(5,'(a70)') TITLE 
  READ(5,*) LREAD,LWRITE,LTEST
  READ(5,*) (LCAL(I),I=1,NPHI)
  READ(5,*) monCell,pRefCell,MPoints
  READ(5,*) SLARGE,SORMAX
  READ(5,*) DENSIT,VISCOS
  READ(5,*) PRANL,TREF,BETA
  READ(5,*) LBUOY,GRAVX,GRAVY,GRAVZ,BOUSSINESQ
  READ(5,*) roughWall,EROUGH,ZZERO
  READ(5,*) FACNAP,FACFLX
  READ(5,*) LTRANSIENT,BDF,BTIME,CN
  READ(5,*) LEVM,LASM,LLES
  READ(5,*) LSGDH,LGGDH,LAFM
  READ(5,*) TurbModel
  READ(5,*) UIN,VIN,WIN,TEIN,EDIN,TIN,VARTIN,CONIN
  READ(5,*) convective_scheme
  READ(5,*) limiter
  READ(5,*) (GDS(I),I=1,NPHI)
  READ(5,*) (URF(I),I=1,NPHI)
  READ(5,*) (SOR(I),I=1,NPHI)
  READ(5,*) (NSW(I),I=1,NPHI)
  READ(5,*) NUMSTEP,TIMESTEP,NZAPIS,MAXIT
  READ(5,*) lstsq, lstsq_qr, lstsq_dm, gauss
  READ(5,*) NPCOR, NIGRAD
  READ(5,*) SIMPLE,PISO,PIMPLE,ncorr
  READ(5,*) const_mflux
  READ(5,*) CoNumFix, CoNumFixValue
!.END: READ INPUT FILE.............................................!
  CLOSE (5)

!.Create an input file reading log:
  WRITE(6,'(a)') '  Input file log: '
  WRITE(6,'(a)') '---cut here-----------------------------------------------------------------------------'
  WRITE(6,'(a70)') TITLE
  WRITE(6,'(3(L1,1x),5x,a)') LREAD,LWRITE,LTEST,'READ3,WRIT3,LTEST'
  WRITE(6,'(10(L1,1x),5x,a)') (LCAL(I),I=1,NPHI),'(LCAL(I),I=1,NPHI),IP=4,ITE=5,IED=6,IEN=7,IVIS=8,IVART=9,ICON=10'
  WRITE(6,'(3(i3,1x),5x,a)') monCell,pRefCell,MPoints,'monCell,pRefCell,MPoints'
  WRITE(6,'(2(es11.4,1x),5x,a)') SLARGE,SORMAX,'SLARGE,SORMAX'
  WRITE(6,'(2(es11.4,1x),a)') DENSIT,VISCOS,'DENSIT,VISCOS'
  WRITE(6,'(3(es11.4,1x),a)') PRANL,TREF,BETA,'PRANL,TREF,BETA'
  WRITE(6,'(L1,1x,3f6.2,1x,i1,1x,a)') LBUOY,GRAVX,GRAVY,GRAVZ,BOUSSINESQ,'LBUOY,GRAVX,GRAVY,GRAVZ,BOUSSINESQ'
  WRITE(6,'(L1,1x,f5.2,1x,es11.4,1x,a)') roughWall,EROUGH,ZZERO,'roughWall,EROUGH,ZZERO'
  WRITE(6,'(2(f4.2,1x),a)') FACNAP,FACFLX,'FACNAP,FACFLX'
  WRITE(6,'(L1,1x,L1,1x,f4.2,1x,L1,1x,a)') LTRANSIENT,BDF,BTIME,CN,'LTRANSIENT,BDF,BTIME,CN'
  WRITE(6,'(3(L1,1x),a)') LEVM,LASM,LLES,'LEVM,LASM,LLES'
  WRITE(6,'(3(L1,1x),a)') LSGDH,LGGDH,LAFM,'LSGDH,LGGDH,LAFM'
  WRITE(6,'(i2,1x,a)') TurbModel, 'Turbulence Model'
  WRITE(6,'(8(es11.4,1x),a)') UIN,VIN,WIN,TEIN,EDIN,TIN,VARTIN,CONIN,'UIN,VIN,WIN,TEIN,EDIN,TIN,VARTIN,CONIN'
  WRITE(6,'(a,a)') convective_scheme, 'Convective scheme'
  WRITE(6,'(a,1x,a)') limiter, 'Gradient limiter'
  WRITE(6,'(10(f4.2,1x),a)') (GDS(I),I=1,NPHI),'(GDS(I),I=1,NPHI)'
  WRITE(6,'(10(f4.2,1x),a)') (URF(I),I=1,NPHI),'(URF(I),I=1,NPHI)'
  WRITE(6,'(10(es9.2,1x),a)') (SOR(I),I=1,NPHI),'(SOR(I),I=1,NPHI)'
  WRITE(6,'(10(i3,1x),a)') (NSW(I),I=1,NPHI),'(NSW(I),I=1,NPHI)'
  WRITE(6,'(i5,1x,es9.2,1x,i5,1x,i4,1x,a)') NUMSTEP,TIMESTEP,NZAPIS,MAXIT,'NUMSTEP,TIMESTEP,NZAPIS,MAXIT'
  WRITE(6,'(4(L1,1x),a)') lstsq, lstsq_qr, lstsq_dm, gauss,'lstsq, lstsq_qr, lstsq_dm, gauss'
  WRITE(6,'(i1,1x,i1,1x,a)') NPCOR, NIGRAD,'NPCOR, NIGRAD'
  WRITE(6,'(3(L1,1x),i1,1x,a)') SIMPLE,PISO,PIMPLE,ncorr,'SIMPLE,PISO,PIMPLE,ncorr'
  WRITE(6,'(1(L1,1x),5x,a)') const_mflux,'const_mflux'
  WRITE(6,'(L1,es11.4,5x,a)') CoNumFix, CoNumFixValue,'CoNumFix, CoNumFixValue'
  WRITE(6,'(a)') '---cut here-----------------------------------------------------------------------------'

  !
  ! Turbulent flow computation condition:
  !
  lturb = levm.or.lasm.or.lles

  !
  ! Convective scheme:
  !
  if(adjustl(convective_scheme) == 'central') then
    lcds = .true.
  elseif(adjustl(convective_scheme) == 'cds-corrected') then
    lcdsc = .true.
  elseif(adjustl(convective_scheme) == 'linear') then
    lluds = .true.
  elseif(adjustl(convective_scheme) == 'smart') then
    lsmart = .true.
  elseif(adjustl(convective_scheme) == 'avl-smart') then
    lavl = .true.
  elseif(adjustl(convective_scheme) == 'muscl') then
    lmuscl = .true.
  elseif(adjustl(convective_scheme) == 'umist') then
    lumist = .true.
  elseif(adjustl(convective_scheme) == 'koren') then
    lkoren = .true.
  elseif(adjustl(convective_scheme) == 'charm') then
    lcharm = .true.
  elseif(adjustl(convective_scheme) == 'ospre') then
    lospre = .true.
  elseif(adjustl(convective_scheme) == 'central-f') then
    lcds_flnt = .true.
  elseif(adjustl(convective_scheme) == 'linear-f') then
    l2nd_flnt = .true.
  elseif(adjustl(convective_scheme) == 'limited-linear') then
    l2ndlim_flnt = .true.
  elseif(adjustl(convective_scheme) == 'muscl-f') then
    lmuscl_flnt = .true.
  else
    write(*,'(a)') '  Convective scheme not chosen, assigning default muscl scheme'
    convective_scheme = 'muscl'
  endif
  
  ! Set value for flux_limiter logical
  if(lluds.or.lsmart.or.lavl.or.lmuscl.or.lumist.or.lkoren.or.lcharm.or.lospre) then
    flux_limiter = .true.
  else
    flux_limiter = .false.
  endif

  write(*,'(a)') ' '
  write(*,'(2a)') '  Convective scheme: ', adjustl(convective_scheme)
  write(*,'(a)') ' '

  !
  ! Gradient limiter:
  !
  if(adjustl(limiter) == 'Barth-Jespersen') then
    write(*,*) ' Gradient limiter: Barth-Jespersen'
  elseif(adjustl(limiter) == 'Venkatakrishnan') then
    write(*,*) ' Gradient limiter: Venkatakrishnan'
  elseif(adjustl(limiter) == 'MVenkatakrishnan') then
    write(*,*) ' Gradient limiter: Wang modified Venkatakrishnan'
  elseif(adjustl(limiter) == 'no-limit') then
    write(*,*) ' Gradient limiter: no-limit'
  else
    write(*,*) ' Gradient limiter type not chosen, assigning default Venkatakrishnan limiter'
    limiter = 'Venkatakrishnan'
  endif

  !
  ! Open files for data at monitoring points 
  !
  if(ltransient) then
    open(unit=89,file=trim(out_folder_path)//'/transient_monitoring_points')
    rewind 89
    do imon=1,mpoints
      write(trpn,'(i2)') imon
      open(91+imon,file=trim(out_folder_path)//"/transient_monitor_point_"//trpn, access='append')
      if(.not.lread) rewind(91+imon)
    end do
  end if

end subroutine