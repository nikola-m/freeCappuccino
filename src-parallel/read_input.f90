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

  include 'mpif.h'

  integer :: i,imon
  character(len=2) :: trpn
  character(len=25) :: convective_scheme
!
!***********************************************************************
!

  ! Root processor opens files
  if (myid .eq. 0) then

  open(unit=5,file=input_file)
  rewind 5

  read(5,'(a70)') title 
  read(5,*) lread,lwrite,ltest
  read(5,*) (lcal(i),i=1,nphi)
  read(5,*) monCell,pRefCell,MPoints
  read(5,*) slarge,sormax
  read(5,*) densit,viscos
  read(5,*) pranl,tref,beta
  read(5,*) lbuoy,gravx,gravy,gravz,boussinesq
  read(5,*) roughWall,EROUGH,ZZERO
  read(5,*) facnap,facflx
  read(5,*) ltransient,bdf,btime,cn
  read(5,*) levm,lasm,lles,ldes
  read(5,*) lsgdh,lggdh,lafm
  read(5,*) TurbModel
  read(5,*) uin,vin,win,tein,edin,tin,vartin,conin
  read(5,*) convective_scheme
  read(5,*) limiter
  read(5,*) (gds(i),i=1,nphi)
  read(5,*) (urf(i),i=1,nphi)
  read(5,*) (sor(i),i=1,nphi)
  read(5,*) (nsw(i),i=1,nphi)
  read(5,*) numstep,timestep,nzapis,maxit
  read(5,*) lstsq, lstsq_qr, lstsq_dm, gauss
  read(5,*) npcor, nigrad
  read(5,*) simple,piso,pimple,ncorr
  read(5,*) const_mflux
  read(5,*) CoNumFix, CoNumFixValue

  close (5)

!.Create an input file reading log:
  write(6,'(a)') '  Input file log: '
  write(6,'(a)') '---cut here-----------------------------------------------------------------------------'
  write(6,'(a70)') title
  write(6,'(3(l1,1x),5x,a)') lread,lwrite,ltest,'read3,writ3,ltest'
  write(6,'(10(l1,1x),5x,a)') (lcal(i),i=1,nphi),'(lcal(i),i=1,nphi),ip=4,ite=5,ied=6,ien=7,ivis=8,ivart=9,icon=10'
  write(6,'(3(i3,1x),5x,a)') monCell,pRefCell,MPoints,'monCell,pRefCell,MPoints'
  write(6,'(2(es11.4,1x),5x,a)') slarge,sormax,'slarge,sormax'
  write(6,'(2(es11.4,1x),a)') densit,viscos,'densit,viscos'
  write(6,'(3(es11.4,1x),a)') pranl,tref,beta,'pranl,tref,beta'
  write(6,'(l1,1x,3f6.2,1x,l1,1x,a)') lbuoy,gravx,gravy,gravz,boussinesq,'lbuoy,gravx,gravy,gravz,boussinesq'
  write(6,'(L1,1x,f5.2,1x,es11.4,1x,a)') roughWall,erough,zzero,'roughWall,erough,zzero'
  write(6,'(2(f4.2,1x),a)') facnap,facflx,'facnap,facflx'
  write(6,'(l1,1x,l1,1x,f4.2,1x,l1,1x,a)') ltransient,bdf,btime,cn,'ltransient,bdf,btime,cn'
  write(6,'(4(l1,1x),a)') levm,lasm,lles,ldes,'levm,lasm,lles,ldes'
  write(6,'(3(l1,1x),a)') lsgdh,lggdh,lafm,'lsgdh,lggdh,lafm'
  write(6,'(i2,1x,a)') TurbModel, 'Turbulence Model'
  write(6,'(8(es11.4,1x),a)') uin,vin,win,tein,edin,tin,vartin,conin,'uin,vin,win,tein,edin,tin,vartin,conin'
  write(6,'(a,a)') convective_scheme, 'Convective scheme'
  write(6,'(a,1x,a)') limiter, 'Gradient limiter'
  write(6,'(10(f4.2,1x),a)') (gds(i),i=1,nphi),'(gds(i),i=1,nphi)'
  write(6,'(10(f4.2,1x),a)') (urf(i),i=1,nphi),'(urf(i),i=1,nphi)'
  write(6,'(10(es9.2,1x),a)') (sor(i),i=1,nphi),'(sor(i),i=1,nphi)'
  write(6,'(10(i3,1x),a)') (nsw(i),i=1,nphi),'(nsw(i),i=1,nphi)'
  write(6,'(i5,1x,es9.2,1x,i5,1x,i4,1x,a)') numstep,timestep,nzapis,maxit,'numstep,timestep,nzapis,maxit'
  write(6,'(4(L1,1x),a)') lstsq, lstsq_qr, lstsq_dm, gauss,'lstsq, lstsq_qr, lstsq_dm, gauss'
  write(6,'(i1,1x,i1,1x,a)') npcor, nigrad,'npcor, nigrad'
  write(6,'(3(l1,1x),i1,1x,a)') simple,piso,pimple,ncorr,'simple,piso,pimple,ncorr'
  write(6,'(1(L1,1x),5x,a)') const_mflux,'const_mflux'
  write(6,'(L1,es11.4,5x,a)') CoNumFix, CoNumFixValue,'CoNumFix, CoNumFixValue'
  write(6,'(a)') '---cut here-----------------------------------------------------------------------------'

  endif


! Broadcast input data to other processes
  call MPI_BCAST(title,70,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lread,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lwrite,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(ltest,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lcal,NPHI,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(monCell,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(pRefCell,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(MPOINTS,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  ! Treba naci kom procesu pripadaju monitoring tacke - pogledaj getpidlm.f kod Sase.

  call MPI_BCAST(slarge,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(sormax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(densit,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(viscos,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(pranl,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(tref,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(beta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lbuoy,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gravx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gravy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gravz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(boussinesq,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(roughWall,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(erough,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(zzero,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(facnap,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(facflx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(ltransient,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(btime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(bdf,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(cn,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(levm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lasm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lles,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(ldes,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(lsgdh,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lggdh,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lafm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(turbmodel,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR) 

  call MPI_BCAST(uin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(vin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(win,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(tein,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(edin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(tin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(vartin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(conin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)


  call MPI_BCAST(convective_scheme,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(limiter,20,MPI_CHARACTER,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(GDS(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(URF(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(SOR(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(NSW(1),NPHI,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(numstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(timestep,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(nzapis,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(maxit,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)


  call MPI_BCAST(lstsq,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lstsq_qr,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(lstsq_dm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(gauss,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(npcor,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(nigrad,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(simple,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(piso,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(pimple,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(ncorr,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(const_mflux,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)

  call MPI_BCAST(conumfix,1,MPI_LOGICAL,0,MPI_COMM_WORLD,IERR)
  call MPI_BCAST(conumfixvalue,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERR)

  ! Izgubi mu se pojam koji je process rank pa moram ovo da pozovem:
  call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr  ) 

if (myid .eq. 0) then
  write(6,*)' '
  write(6,*)'  Finished reading and broadcasting input data.'
  write(6,*)' '
endif


  !
  ! Turbulent flow computation condition:
  !
  lturb = levm.or.lasm.or.lles.or.ldes

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
  elseif(adjustl(convective_scheme) == 'muscl-f') then
    lmuscl_flnt = .true.
  else
    if (myid .eq. 0) then
      write(*,'(a)') '  Convective scheme not chosen, assigning default muscl scheme'
    endif
    convective_scheme = 'muscl'
  endif
  
  ! Set value for flux_limiter logical
  if(lluds.or.lsmart.or.lavl.or.lmuscl.or.lumist.or.lkoren.or.lcharm.or.lospre) then
    flux_limiter = .true.
  else
    flux_limiter = .false.
  endif

  if (myid .eq. 0) then
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

    elseif(adjustl(limiter) == 'mVenkatakrishnan') then

      write(*,*) ' Gradient limiter: Wang modified Venkatakrishnan'

    else!if(adjustl(limiter) == 'no-limit') then

      write(*,*) ' Gradient limiter: no-limit'
      
    endif

    write(*,'(a)') ' '

  !
  ! Time stepping algorithm:
  !
    if( bdf ) then

      if (btime < 1.) then
        write(*,*) ' Time stepping method: Euler Implicit'
      else
        write(*,*) ' Time stepping method: Three Level Implicit Time Integration (BDF2)'
      endif

    elseif( cn ) then

         write(*,*) ' Time stepping method: Crank-Nicolson'   

    endif

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

  !
  ! Pressure reference cell
  !

  iPrefProcess = 0
  pRefCell = 1
  
end subroutine