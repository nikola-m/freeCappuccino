!***********************************************************************
!
subroutine init
!
!***********************************************************************
!     Contents:
!
! 0)  Print code logo and timestamp in monitor file
! 1)  Open & Read Input File
! 2)  Open & Read Grid File & Allocate Arrays
! 3)  Index arrays of matrix elements stored in CSR format
! 4)  Set Index Arrays For Cell Looping, Set Monitoring Point And Pressure Reference Point
! 5)  Various initialisations
!     5.1)  Field Initialisation
! 6)  Read Restart File And Set Field Values
! 7)  Initial Gradient Calculation
! 8)  Initialization of parameters for boundary adjecent cells
! 9)  Calculate distance to the nearest wall.
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use variables
  use title_mod
  use gradients
  use sparse_matrix, only: create_CSR_matrix_from_mesh_data
  use k_epsilon_std, only: te,ed,dTEdxi,dEDdxi,allocate_k_epsilon_std
  use temperature, only: t,utt,vtt,wtt,pranl
  use utils, only: timestamp, show_logo, i4vec_print2

  implicit none
!
!***********************************************************************
!

! 
! Local variables 
!
  integer :: i, ijp, inp, iface
  real(dp) :: nxf, nyf, nzf
  real(dp) :: are

!
!***********************************************************************
!

! 
! 0)  Print code logo and timestamp in monitor file
!
      call show_logo

!
! 1)  Open & Read Input File
!

!.....OPEN & READ INPUT FILE...................................................
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
      READ(5,*) LTRANSIENT,BTIME
      READ(5,*) LEVM,LASM,LLES
      READ(5,*) LSGDH,LGGDH,LAFM
      READ(5,*) TurbModel
      READ(5,*) UIN,VIN,WIN,TEIN,EDIN,TIN,VARTIN,CONIN
      READ(5,*) LCDS,LLUDS,LSMART,LAVL,LMUSCL,LUMIST,LGAMMA
      READ(5,*) (GDS(I),I=1,NPHI)
      READ(5,*) (URF(I),I=1,NPHI)
      READ(5,*) (SOR(I),I=1,NPHI)
      READ(5,*) (NSW(I),I=1,NPHI)
      READ(5,*) NUMSTEP,TIMESTEP,NZAPIS,MAXIT
      READ(5,*) lstsq, lstsq_qr, lstsq_dm, gauss
      READ(5,*) NPCOR, NIGRAD
      READ(5,*) BDF,CN
      READ(5,*) SIMPLE,PISO,PIMPLE,ncorr
      READ(5,*) const_mflux
      READ(5,*) CoNumFix, CoNumFixValue
!.....END: READ INPUT FILE.............................................!
      CLOSE (5)

!.....Create an input file reading log:
      WRITE(6,'(a)') ' Input file log: '
      WRITE(6,'(a)') '--------------------------------------------------------------------------------'
      WRITE(6,'(a70)') TITLE
      WRITE(6,'(3(L1,1x),5x,a)') LREAD,LWRITE,LTEST,'READ3,WRIT3,LTEST'
      WRITE(6,'(10(L1,1x),5x,a)') (LCAL(I),I=1,NPHI),'(LCAL(I),I=1,NPHI),IP=4,ITE=5,IED=6,IEN=7,IVIS=8,IVART=9,ICON=10'
      WRITE(6,'(3(i3,1x),5x,a)') monCell,pRefCell,MPoints,'monCell,pRefCell,MPoints'
      WRITE(6,'(2(es11.4,1x),5x,a)') SLARGE,SORMAX,'SLARGE,SORMAX'
      WRITE(6,'(2(es11.4,1x),a)') DENSIT,VISCOS,'DENSIT,VISCOS'
      WRITE(6,'(3(es11.4,1x),a)') PRANL,TREF,BETA,'PRANL,TREF,BETA'
      WRITE(6,'(L1,1x,3f5.2,1x,i1,1x,a)') LBUOY,GRAVX,GRAVY,GRAVZ,BOUSSINESQ,'LBUOY,GRAVX,GRAVY,GRAVZ,BOUSSINESQ'
      WRITE(6,'(L1,1x,f5.2,1x,es11.4,1x,a)') roughWall,EROUGH,ZZERO,'roughWall,EROUGH,ZZERO'
      WRITE(6,'(2(f4.2,1x),a)') FACNAP,FACFLX,'FACNAP,FACFLX'
      WRITE(6,'(L1,1x,f4.2,1x,a)') LTRANSIENT,BTIME,'LTRANSIENT,BTIME'
      WRITE(6,'(3(L1,1x),a)') LEVM,LASM,LLES,'LEVM,LASM,LLES'
      WRITE(6,'(3(L1,1x),a)') LSGDH,LGGDH,LAFM,'LSGDH,LGGDH,LAFM'
      WRITE(6,'(i2,1x,a)') TurbModel, 'TurbModel'
      WRITE(6,'(8(es11.4,1x),a)') UIN,VIN,WIN,TEIN,EDIN,TIN,VARTIN,CONIN,'UIN,VIN,WIN,TEIN,EDIN,TIN,VARTIN,CONIN'
      WRITE(6,'(7(L1,1x),a)') LCDS,LLUDS,LSMART,LAVL,LMUSCL,LUMIST,LGAMMA,'LCDS,LLUDS,LQUDS,LSMART,LAVL,LMUSCL,LUMIST,LGAMMA'
      WRITE(6,'(10(f4.2,1x),a)') (GDS(I),I=1,NPHI),'(GDS(I),I=1,NPHI), MUSCL velocity, CDS other'
      WRITE(6,'(10(f4.2,1x),a)') (URF(I),I=1,NPHI),'(URF(I),I=1,NPHI)'
      WRITE(6,'(10(es9.2,1x),a)') (SOR(I),I=1,NPHI),'(SOR(I),I=1,NPHI)'
      WRITE(6,'(10(i3,1x),a)') (NSW(I),I=1,NPHI),'(NSW(I),I=1,NPHI)'
      WRITE(6,'(i5,1x,es9.2,1x,i5,1x,i4,1x,a)') NUMSTEP,TIMESTEP,NZAPIS,MAXIT,'NUMSTEP,TIMESTEP,NZAPIS,MAXIT'
      WRITE(6,'(4(L1,1x),a)') lstsq, lstsq_qr, lstsq_dm, gauss,'lstsq, lstsq_qr, lstsq_dm, gauss'
      WRITE(6,'(i1,1x,i1,1x,a)') NPCOR, NIGRAD,'NPCOR, NIGRAD'
      WRITE(6,'(2(L1,1x),1x,a)') BDF,CN,'BDF,CN'
      WRITE(6,'(3(L1,1x),i1,1x,a)') SIMPLE,PISO,PIMPLE,ncorr,'SIMPLE,PISO,PIMPLE,ncorr'
      WRITE(6,'(1(L1,1x),5x,a)') const_mflux,'const_mflux'
      WRITE(6,'(L1,es11.4,5x,a)') CoNumFix, CoNumFixValue,'CoNumFix, CoNumFixValue'
      WRITE(6,'(a)') '--------------------------------------------------------------------------------'
      WRITE(6,'(a)') ' '


!
! 2)  Open & Read mesh file, calculate mesh geometrical quantities, allocate arrays
!
  call mesh_geometry

  call set_parameters  

  call allocate_arrays

  call allocate_gradients

  select case (TurbModel)
    case (1)
      call allocate_k_epsilon_std
    case default
  end select  

!
! 3)  Index arrays of matrix elements stored in CSR format
!
  call create_CSR_matrix_from_mesh_data

!
! 4)  Set Coefficient values for Turbulence models
!

  ! Turbulent flow computation
  lturb = levm.or.lles


  ! Reciprocal values of underrrelaxation factors
  do i=1,nphi
    urfr(i)=1.0_dp / urf(i)
    urfm(i)=1.0_dp - urf(i)
  enddo




!
! 5)  Various initialisations
!

! 5.0)  Parameter Initialisation

  ! Initial time
  if(.not.lread) time=0.0d0

  !Set to zero cumulative error in continuity
  cumulativeContErr = 0.0_dp

! Bulk velocity - important const_mflux flow!
  magUbar = uin

! 5.1)  Field Initialisation - Turbulence kinetic energy and dissipation at inlet

  ! IF(Wilcox.or.SST.or.SAS.or.EARSM_WJ.or.EARSM_M) THEN
  !   TEIN=1e-6*UIN*UIN
  !   EDIN=5.*UIN/(zc(intc)-zc(inbc))
  ! ELSEIF(stdkeps.or.Durbin.or.RNG.or.Realizable.or.LowRe_LB) THEN 
  !   TEIN=1.5*(0.05*UIN)**2 ! for Ti=0.05  
  !   EDIN=cmu75*tein**1.5/(zc(intc)-zc(inbc)) 
  ! ENDIF

! 5.2)  Field Initialisation - reading inlet file


! 5.3)  Field Initialisation

! Field initialisation loop over inner cells--------------------------------
  do inp = 1,numCells

! Initialization of field variables from input file:
  u(inp) = xc(inp)+yc(inp)+zc(inp)!uin
  v(inp) = vin
  w(inp) = win
  te(inp) = tein
  ed(inp) = edin

! Channel flow:
! Random number based fluctuation of mean profile            
      ! call init_random_seed()
      ! call random_number(perturb)    
      ! perturb = 0.9+perturb/5. ! Max perturbation is +/- 10% of mean profile
      ! u(inp) = perturb*u(inp)
      ! v(inp) = perturb/100.
      ! w(inp) = perturb/100.

  enddo

  do i=1,numBoundaryFaces
  iface = numInnerFaces + i
  ijp = numCells+i
  u(ijp) = xf(iface)+yf(iface)+zf(iface)
  enddo

  !-------------------------------------------------------    
  ! Field initialisation over inner cells + boundary faces
  !-------------------------------------------------------

  ! Density
  den = densit
  ! Effective viscosity
  vis = viscos
  ! Temperature
  ! t = tin
  ! ! Temperature variance
  ! vart=vartin
  ! ! Concentration
  ! con=conin

  ! Reynolds stress tensor components
  uu = 0.0_dp
  vv = 0.0_dp
  ww = 0.0_dp
  uv = 0.0_dp
  uw = 0.0_dp
  vw = 0.0_dp

  ! ! Turbulent heat fluxes
  ! utt = 0.0_dp
  ! vtt = 0.0_dp
  ! wtt = 0.0_dp

  ! Reynolds stress anisotropy
  ! if(earsm_wj.or.earsm_m) bij(:,:)=0.0_dp

  ! Pressure and pressure correction
  p = small
  pp = p



!
! 6)  Read Restart File And Set Field Values
!
  if(lread) then
    call readfiles
    pp = p
  end if



!
! 7)  Initial Gradient Calculation
!
  dUdxi = 0.0_dp
  dVdxi = 0.0_dp
  dWdxi = 0.0_dp
  dPdxi = 0.0_dp
  ! dTEdxi = 0.0_dp
  ! dEDdxi = 0.0_dp


  if (lstsq .or. lstsq_qr .or. lstsq_dm) then
    call create_lsq_gradients_matrix(U,dUdxi)
  endif

  call grad(U,dUdxi)
  ! call grad(V,dVdxi)
  ! call grad(W,dWdxi)

! print*,'gradijenti:'
! do i=1,numCells
!   print*,i,':',u(i),dudxi(1,i)
! enddo




! 8) Calculate distance dnw of wall adjecent cells and distance to the nearest wall of all cell centers.

  ! Loop over wall boundaries to calculate normal distance from cell center dnw.
  do i = 1,nwal
    iface = iWallFacesStart+i
    ijp = owner(iface)

    ! Face area 
    are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

    ! Face normals
    nxf = arx(iface)/are
    nyf = ary(iface)/are
    nzf = arz(iface)/are

    ! We need the minus sign because of the direction of normal vector to boundary face which is positive if it faces out.
    dnw(i) = (xf(iface)-xc(ijp))*nxf + (yf(iface)-yc(ijp))*nyf + (zf(iface)-zc(ijp))*nzf

    ! Cell face area divided by distance to the cell center
    srdw(i) = are/dnw(i)
  enddo

  ! Loop over symmetry boundaries to calculate normal distance from cell center dns.
  do i = 1,nsym
    iface = iSymmetryFacesStart+i
    ijp = owner(iface)

    ! Face area 
    are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

    ! Face normals
    nxf = arx(iface)/are
    nyf = ary(iface)/are
    nzf = arz(iface)/are

    ! We need the minus sign because of the direction of normal vector to boundary face which is positive if it faces out.
    dns(i) = (xf(iface)-xc(ijp))*nxf + (yf(iface)-yc(ijp))*nyf + (zf(iface)-zc(ijp))*nzf

    ! Cell face area divided by distance to the cell center
    srds(i) = are/dns(i)

  enddo

  ! srdoc

  ! foc


  !
  ! 9) Distance to the nearest wall needed for some turbulence models
  !

      ! Source term
!      do k=2,nkm; do i=2,nim; do j=2,njm
!      inp=lk(k)+li(i)+j
!          ! Wall distance Poisson eq. source :
!          su(inp) = Vol(inp) 
!      enddo; enddo; enddo  





      ! ! Initialize solution and set fixedValue boundaries
      ! pp=0. 
      ! do k=1,nk; do i=1,ni; do j=1,nj
      ! inp=lk(k)+li(i)+j
      ! intc = lk(nk)+li(i)+j
      ! inbc = lk(1)+li(i)+j
      !     ! 1. If channel - upper and lower boundaries are Wall:
      !     pp(inp) = min(zc(intc)-zc(inp),zc(inp)-zc(inbc))
      !     ! 2. If Boundary layer - only lower boundary is Wall:
      !     !pp(inp) = zc(inp)-zc(inbc)
      ! enddo; enddo; enddo   

      ! wallDistance = pp
      ! pp=p





      ! Initial internal field
      !do k=2,nkm; do i=2,nim; do j=2,njm
      !inp=lk(k)+li(i)+j
      !    pp(inp) = 0.
      !enddo; enddo; enddo 
 
!      sv = 1.0_dp                ! Unit coefficient array
!      call fvm_laplacian(sv,pp) ! Laplacian operator and BCs

!      call cgstab_sip(pp,ip) 

      ! Gradient of solution field stored in pp (gradient stored in dPdxi) :
!      call grad_lsq_qr(pp,dPdxi,2,d)

      ! Wall distance computation from Poisson eq. solution stored in pp:
!      wallDistance = -sqrt(  dPdxi(1,:)*dPdxi(1,:)+dPdxi(2,:)*dPdxi(2,:)+dPdxi(3,:)*dPdxi(3,:)  ) + &
!                      sqrt(  dPdxi(1,:)*dPdxi(1,:)+dPdxi(2,:)*dPdxi(2,:)+dPdxi(3,:)*dPdxi(3,:) + 2.*pp  )

 
!      su = 0.0_dp; sv = 0.0_dp ! Clear arrays

!      call plot_3D_field_vtk (88, trim(out_folder_path)//'/wallDistance_scalar_field', 'scalar', &
!                             'vtk', 'WDIS_field', 'wall-distance ',                      &
!                              NI, NJ, NK, 1, 1, 1, Xc, Yc, Zc, wallDistance, 0.0, 0.0)

!       Open(Unit=87,File=Trim(Out_Folder_Path)//'/wallDistance.plt') 
!       Rewind 87
!       Write(87,*) 'Title     = " "'
!       Write(87,*) 'Variables = "X"'
!       Write(87,*) '"Y"'
!       Write(87,*) '"Z"'
!       Write(87,*) '"Wdist"'
!       Write(87,*) 'Zone T=" "'
!       Write(87,*) 'I=',Ni, ' ,J=',Nj, ' ,K=',Nk,', F=Point'
!       Do k=1,nk; do j=1,nj; do i=1,ni
!       Inp=Lk(K)+Li(I)+J
!       Write(87,*) Xc(Inp),Yc(Inp),Zc(Inp),wallDistance(Inp)
!       Enddo; Enddo; Enddo 
!       Close(87)        

end subroutine
