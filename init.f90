!***********************************************************************
!
subroutine init
!
!***********************************************************************
!     Contents:
!
! 0)  Print code logo and timestamp in monitor file
! 1)  Open & Read Input File
! 2)  Open & Read Grid File & Allocating Arrays
! 3)  Index arrays of matrix elements stored in CSR format
! 4)  Set Index Arrays For Cell Looping, Set Monitoring Point And Pressure Reference Point
! 5)  Set Initial timestepping control values
! 6)  Various initialisations
!     6.1)  Field Initialisation
! 7)  Read Restart File And Set Field Values
! 8)  Initial Gradient Calculation
! 9) Calculate distance to the nearest wall.
!
!***********************************************************************
  use types
  use parameters
  use indexes
  use geometry
  use variables
  use title_mod
  use gradients
  use sparse_matrix, only: create_CSR_matrix_from_mesh_data
  use k_epsilon_std, only: te,ed,dTEdxi,dEDdxi,allocate_k_epsilon_std
  use temperature, only: t,utt,vtt,wtt,pranl
  use utils

  implicit none
!
!***********************************************************************
!

! 
! Local variables 
!
  integer :: i, inp
  real(dp) :: perturb

!
!***********************************************************************
!

! 
! 0)  Print code logo and timestamp in monitor file
!
      call show_logo
      call timestamp ()


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
      READ(5,*) DENSIT,VISCOS,RAY
      READ(5,*) PRANL,TREF,BETA
      READ(5,*) LBUOY,GRAVX,GRAVY,GRAVZ,BOUSSINESQ
      READ(5,*) roughWall,EROUGH,ZZERO
      READ(5,*) PHIT,SKSI,ETA,RCOST,FACNAP,FACFLX
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
      READ(5,*) const_mflux,gradPcmf
      READ(5,*) CoNumFix, CoNumFixValue
!.....END: READ INPUT FILE.............................................!
      CLOSE (5)

!.....Create an input file reading log:
      WRITE(66,'(a)') ' Input file log: '
      WRITE(66,'(a)') '--------------------------------------------------------------------------------'
      WRITE(66,'(a70)') TITLE
      WRITE(66,'(3(L1,1x),5x,a)') LREAD,LWRITE,LTEST,'READ3,WRIT3,LTEST'
      WRITE(66,'(10(L1,1x),5x,a)') (LCAL(I),I=1,NPHI),'(ICAL(I),I=1,NPHI),IEN=7!!!,IVIS=8,IVART=9;ICON=10!'
      WRITE(66,'(3(i3,1x),5x,a)') monCell,pRefCell,MPoints,'monCell,pRefCell,MPoints'
      WRITE(66,'(2(es11.4,1x),5x,a)') SLARGE,SORMAX,'SLARGE,SORMAX'
      WRITE(66,'(3(es11.4,1x),a)') DENSIT,VISCOS,RAY,'DENSIT,VISCOS,RAY'
      WRITE(66,'(3(es11.4,1x),a)') PRANL,TREF,BETA,'PRANL,TREF,BETA'
      WRITE(66,'(L1,1x,3f5.2,1x,i1,1x,a)') LBUOY,GRAVX,GRAVY,GRAVZ,BOUSSINESQ,'LBUOY,GRAVX,GRAVY,GRAVZ,BOUSSINESQ'
      WRITE(66,'(L1,1x,f5.2,1x,es11.4,1x,a)') roughWall,EROUGH,ZZERO,'roughWall,EROUGH,ZZERO'
      WRITE(66,'(6(f4.2,1x),a)') PHIT,SKSI,ETA,RCOST,FACNAP,FACFLX,'PHIT,SKSI,ETA,RCOST,FACNAP,FACFLX'
      WRITE(66,'(L1,1x,f4.2,1x,a)') LTRANSIENT,BTIME,'LTRANSIENT,BTIME'
      WRITE(66,'(3(L1,1x),a)') LEVM,LASM,LLES,'LEVM,LASM,LLES'
      WRITE(66,'(3(L1,1x),a)') LSGDH,LGGDH,LAFM,'LSGDH,LGGDH,LAFM'
      WRITE(66,'(i2,1x,a)') TurbModel, 'TurbModel'
      WRITE(66,'(8(es11.4,1x),a)') UIN,VIN,WIN,TEIN,EDIN,TIN,VARTIN,CONIN,'UIN,VIN,WIN,TEIN,EDIN,TIN,VARTIN,CONIN'
      WRITE(66,'(7(i2,1x),a)') LCDS,LLUDS,LSMART,LAVL,LMUSCL,LUMIST,LGAMMA,'LCDS,LLUDS,LQUDS,LSMART,LAVL,LMUSCL,LUMIST,LGAMMA'
      WRITE(66,'(10(f4.2,1x),a)') (GDS(I),I=1,NPHI),'(GDS(I),I=1,NPHI), MUSCL velocity, CDS other'
      WRITE(66,'(10(f4.2,1x),a)') (URF(I),I=1,NPHI),'(URF(I),I=1,NPHI)'
      WRITE(66,'(10(es9.2,1x),a)') (SOR(I),I=1,NPHI),'(SOR(I),I=1,NPHI)'
      WRITE(66,'(10(i3,1x),a)') (NSW(I),I=1,NPHI),'(NSW(I),I=1,NPHI)'
      WRITE(66,'(i5,1x,es9.2,1x,i5,1x,i4,1x,a)') NUMSTEP,TIMESTEP,NZAPIS,MAXIT,'NUMSTEP,TIMESTEP,NZAPIS,MAXIT'
      WRITE(66,'(4(L1,1x),a)') lstsq, lstsq_qr, lstsq_dm, gauss,'lstsq, lstsq_qr, lstsq_dm, gauss'
      WRITE(66,'(i1,1x,i1,1x,a)') NPCOR, NIGRAD,'NPCOR, NIGRAD'
      WRITE(66,'(2(L1,1x),1x,a)') BDF,CN,'BDF,CN'
      WRITE(66,'(3(L1,1x),i1,1x,a)') SIMPLE,PISO,PIMPLE,ncorr,'SIMPLE,PISO,PIMPLE,ncorr'
      WRITE(66,'(1(L1,1x),es11.4,5x,a)') const_mflux,gradPcmf,'periodic_boundary, const_mflux, gradPcmf'
      WRITE(66,'(L1,es11.4,5x,a)') CoNumFix, CoNumFixValue,'CoNumFix, CoNumFixValue'
      WRITE(66,'(a)') '--------------------------------------------------------------------------------'
      WRITE(66,'(a)') ' '


!
! 2)  Open & Read Grid File & Allocating Arrays
!

      ! Initialize number of boudnary faces for each boundary type 
      ninl = 0 
      nout = 0
      nsym = 0
      nwal = 0
      npru = 0
      noc = 0

      
      open(unit=4,file=grid_file,form='unformatted')
      rewind 4

      read(4) &
            numNodes,numCells,numInnerFaces, &
            ninl,nout,nsym,npru,nwal,noc, &
            nwali,nwala,nwalf


!-----------------------------------------------------------------------
! STOP READING FILE FOR A SECOND AND DO SOME USEFUL WORK
!-----------------------------------------------------------------------
      call set_parameters                                              
      call allocate_arrays
      call allocate_gradients
      select case (TurbModel)
        case (1)
          call allocate_k_epsilon_std
        case default
      end select                                           
!-----------------------------------------------------------------------
 
      if(ninl.gt.0) then     
      read(4) &
              (xni(i),i=1,ninl),(yni(i),i=1,ninl),(zni(i),i=1,ninl), &
              (xfi(i),i=1,ninl),(yfi(i),i=1,ninl),(zfi(i),i=1,ninl)
      endif

      if(nout.gt.0) then 
      read(4) &
              (xno(i),i=1,nout),(yno(i),i=1,nout),(zno(i),i=1,nout), &
              (xfo(i),i=1,nout),(yfo(i),i=1,nout),(zfo(i),i=1,nout)
      endif

      if(nwal.gt.0) then 
      read(4) &
            (srdw(i),i=1,nwal),(dnw(i),i=1,nwal), &
            (xnw(i),i=1,nwal),(ynw(i),i=1,nwal),(znw(i),i=1,nwal), &
            (xfw(i),i=1,nwal),(yfw(i),i=1,nwal),(zfw(i),i=1,nwal)
      endif

      if(nsym.gt.0) then 
      read(4) &
            (srds(i),i=1,nsym),(dns(i),i=1,nsym), &
            (xns(i),i=1,nsym),(yns(i),i=1,nsym),(zns(i),i=1,nsym), &
            (xfs(i),i=1,nsym),(yfs(i),i=1,nsym),(zfs(i),i=1,nsym)
      endif

      if(npru.gt.0) then 
      read(4) &
            (xnpr(i),i=1,npru),(ynpr(i),i=1,npru),(znpr(i),i=1,npru), &
            (xfpr(i),i=1,npru),(yfpr(i),i=1,npru),(zfpr(i),i=1,npru)
      endif

      if(noc.gt.0) then 
      read(4) &
            (ijl(i) ,i=1,noc),(ijr(i)  ,i=1,noc), &
            (srdoc(i), i=1,noc),(foc(i), i=1,noc), &
            (xnoc(i),i=1,noc),(ynoc(i),i=1,noc),(znoc(i),i=1,noc), &
            (xfoc(i),i=1,noc),(yfoc(i),i=1,noc),(zfoc(i),i=1,noc)
      endif

      read(4) &
!
!           Node data: coordinates of vertices
!
            (x(i) ,i=1,numNodes),(y(i) ,i=1,numNodes),(z(i) ,i=1,numNodes), &
!
!           Cell data: cell centers and volumes
!
            (xc(i) ,i=1,numCells),(yc(i) ,i=1,numCells),(zc(i) ,i=1,numCells), &
            (vol(i),i=1,numCells), &
!
!           The owner and neighbour index arrays, interpolation factors for inner faces
!
            (owner(i),i=1,numFacesTotal),(neighbour(i),i=1,numInnerFaces), &
            (facint(i),i=1,numInnerFaces), &
!
!           Face normal vector - its components are face area projections
!
            (arx(i),i=1,numInnerFaces),(ary(i),i=1,numInnerFaces),(arz(i),i=1,numInnerFaces), &
!
!           Face centers for inner faces
!
            (xf(i),i=1,numInnerFaces),(yf(i),i=1,numInnerFaces),(zf(i),i=1,numInnerFaces)

!.....Close mesh file
      close (4)


!
! 3)  Index arrays of matrix elements stored in CSR format
!
  call create_CSR_matrix_from_mesh_data




!
! 4)  Set Coefficient values for Turbulence models
!

      ! Turbulent flow computation
      LTURB=LEVM.OR.LLES

! !=====Coefficients for Sasa's algebraic flux model
!       C1asm = 1.8_dp
!       C2asm = 0.6_dp 
!       C3asm = 0.6_dp

! !=====STANDARD k-epsilon=============================
!       IF (STDKEPS) THEN
!       CMU = 0.09_dp   
!       C1 = 1.44_dp
!       C2 = 1.92_dp
!       C3 = 1.44_dp
! !=====STANDARD k-epsilon=============================
!       ENDIF

! !=====RENORMALIZATION GROUP (RNG) k-epsilon=============================
!       IF (RNG) THEN
!       CMU = 0.0845   ! za RNG model!!!
!       C1 = 1.42      ! za RNG model!!!
!       C2 = 1.68      ! za RNG model!!!
! !=====RENORMALIZATION GROUP (RNG) k-epsilon=============================
!       ENDIF

! !=====REALIZABLE k-epsilon==============================================
!       IF (REALIZABLE) THEN
!       C1 = 1.44      ! za Realizable k-eps model !!!
!       C2 = 1.9       ! za Realizable k-eps model !!!
! !=====END:REALIZABLE k-epsilon==========================================
!       END IF

! !====================================================
! !    Define turbulence model constants here.
! !    k-omega model: Sigma_k=2.0; Sigma_omega=2.0
! !    REFERENCES:
! !    Wilcox1998:
! !    * Wilcox, D. C., "Reassessment of the Scale-Determining Equation for Advanced Turbulence Models," AIAA Journal, Vol. 26, No. 11, 1988, pp. 1299-1310.
! !    * Wilcox, D. C., Turbulence Modeling for CFD, 1st edition, DCW Industries, Inc., La Canada CA, 1993. 
! !    * Wilcox, D. C., Turbulence Modeling for CFD, 2nd edition, DCW Industries, Inc., La Canada CA, 1998. <- This is where Wilcox1998 formulation comes from.
! !    Wilcox2006:
! !    * Wilcox, D. C., "Formulation of the k-omega Turbulence Model Revisited," AIAA Journal, Vol. 46, No. 11, 2008, pp. 2823-2838.
! !    * Wilcox, D. C., Turbulence Modeling for CFD, 3rd edition, DCW Industries, Inc., La Canada CA, 2006. 
! !====================================================
!       IF (Wilcox) THEN
! !.....Wilcox1998 (These are in Fluent 13):
!       ALPHA=13./25.
!       BETTA=0.072
!       BETTAST=0.09
! !.....Wilcox2006:
! !%      ALPHA=13./25.
! !%      BETTA=0.0708
! !%      BETTAST=0.09
! !=====================================================
!       ENDIF
! !====================================================
! !     Define SST, ad SAS-SST model constants.
! !     REFERENCES:
! !     * ANSYS FLUENT THheory Guide p.71
! !     * Menter, F. R., "Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications," AIAA Journal, Vol. 32, No. 8, August 1994, pp. 1598-1605. 
! !     * Menter, F. R., Kuntz, M., and Langtry, R., "Ten Years of Industrial Experience with the SST Turbulence Model," Turbulence, Heat and Mass Transfer 4, ed: K. Hanjalic, Y. Nagano, and M. Tummers, Begell House, Inc., 2003, pp. 625 - 632. 
! !=====================================================
!       IF (SST.OR.SAS) THEN
!       SIGMK1=1./1.176_dp
!       SIGMK2=1.0_dp
!       SIGMOM1=1./2.0_dp
!       SIGMOM2=1./1.168_dp
!       BETAI1=0.075
!       BETAI2=0.0828
!       A1=0.31_dp
! !.....SST-1994 coefficients
! !%      ALPHA1=(BETAI1/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMAOM1)
! !%      ALPHA2=(BETAI2/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMAOM2)
! !.....SST-2003 coefficients. The first is higher than the original constant
! !     definition by approximately 0.43%, and the second is lower by less than 0.08%. 
!       ALPHA1=5./9.
!       ALPHA2=0.44
!       BETTAST=0.09
!       END IF

!       IF (EARSM_WJ) THEN
!       SIGMK1=1.1_dp
!       SIGMK2=1.1_dp
!       SIGMOM1=0.53_dp
!       SIGMOM2=1.0
!       BETAI1=0.0747
!       BETAI2=0.0828
!       A1=0.31_dp
!       ALPHA1=0.518_dp
!       ALPHA2=0.44
!       BETTAST=0.09
!       END IF

!       IF (EARSM_M) THEN
!       SIGMK1=0.5_dp
!       SIGMK2=1.1_dp
!       SIGMOM1=0.5_dp
!       SIGMOM2=0.856_dp
!       BETAI1=0.075
!       BETAI2=0.0828
!       A1=0.31_dp
!       BETTAST=0.09
!       ALPHA1=(BETAI1/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMOM1)
!       ALPHA2=(BETAI2/BETTAST)-CAPPA**2/(SQRT(BETTAST)*SIGMOM2)
!       END IF

!       CMU25=DSQRT(DSQRT(CMU))
!       CMU75=CMU25**3



! ! Prandtl numbers

! !==========================
! !.....PRANDTL NUMBER FOR FLUID
! !.....PRANL FOR WATER, PRANL=7.0
! !.....PRANL FOR AIR  , PRANL=0.7
! !==========================
!       PRM1=1./PRANL
!       PRANT=0.86
!       PRT1=1./PRANT
! !
! !.....RECIPROCAL VALUES OF PRANDTL NUMBERS
!       DO I=1,NPHI
!       PRTINV(I)=1.0_dp
!       END DO

! !=====TURBULENT MODEL SIGMA's===========================================
!       IF(STDKEPS.OR.DURBIN) THEN
! !.....[Standard k-epsilon Coefficient: ]
!       PRTINV(IED)=1./1.3_dp
! !     [Beljaars (Askervein hill paper), 1987 Coefficients: ]
!       !PRTINV(IED)=1./1.85
!       ENDIF

! !=====RENORMALIZATION GROUP (RNG) k-epsilon=============================
!       IF (RNG) THEN
!       PRTINV(ITE)=1./0.7194_dp
!       PRTINV(IED)=1./0.7194_dp
! !=====END:RENORMALIZATION GROUP (RNG) k-epsilon=========================
!       ENDIF

! !=====REALIZABLE k-epsilon==============================================
!       IF (REALIZABLE) THEN
!       PRTINV(ITE) = 1.0_dp
!       PRTINV(IED) = 1./1.20_dp
! !=====END:REALIZABLE k-epsilon==========================================
!       ENDIF

! !=====Wilcox k-omega==============================================
! !.....Wilcox1998:
!       IF (Wilcox) THEN
!       PRTINV(IED)=0.5
!       PRTINV(ITE)=0.5
! !.....Wilcox2006:
! !%      PRTINV(IED)=0.5
! !%     PRTINV(ITE)=0.6
! !=====END:Wilcox k-omega==========================================
!       ENDIF

! !.....Prandtl-Schmidt number for energy equation:
!       PRTINV(IEN)=1./PRANL
! !      IF(LTURB) PRTINV(IEN)=1./PRANT
! !@      PRRAT=PRANL/PRANT
! !@      PFUN=9.24*(PRRAT**0.75-1.0)*(1.0+0.28*EXP(-0.007*PRRAT))
! !@      WRITE(66,*)'pfun: ',pfun

  ! Reciprocal values of underrrelaxation factors
  do i=1,nphi
    urfr(i)=1.0_dp / urf(i)
    urfm(i)=1.0_dp - urf(i)
  enddo







! 5)  Set Initial timestepping control values
      itime=0
      time=-timestep
      !Set to zero cumulative error in continuity
      cumulativeContErr = 0.0_dp



! 6)  Various initialisations

!     6.0)  Field parameter Initialisation

!     # Bulk velocity - estimated from Reynolds number, and written in 'input' file
      magUbar = UIN !!!<----veoma bitno za const_mflux flow!

! !     6.1)  Field Initialisation - Turbulence kinetic energy and dissipation at inlet

!       inbc = li(1)+1
!       intc = lk(nk)+li(1)+1
!       IF(Wilcox.or.SST.or.SAS.or.EARSM_WJ.or.EARSM_M) THEN
!         TEIN=1e-6*UIN*UIN
!         ! 0.056894 = 0.08948-0.028 -> visina kanala!
!         EDIN=5.*UIN/(zc(intc)-zc(inbc))
!       ELSEIF(stdkeps.or.Durbin.or.RNG.or.Realizable.or.LowRe_LB) THEN 
!         TEIN=1.5*(0.05*UIN)**2 ! for Ti=0.05  
!         EDIN=cmu75*tein**1.5/(zc(intc)-zc(inbc)) 
!       ENDIF

!     6.2)  Field Initialisation - reading inlet file


!     6.3)  Field Initialisation

!-----Field initialisation loop over inner cells--------------------------------
      do inp=1,numCells

!.....Initialization of field variables from input file:
      u(inp)=uin
      v(inp)=vin
      w(inp)=win
      te(inp)=tein
      ed(inp)=edin

!.....Channel flow:
!.....Random number based fluctuation of mean profile            
      call init_random_seed()
      call random_number(perturb)    
      perturb = 0.9+perturb/5. ! Max perturbation is +/- 10% of mean profile
      u(inp) = perturb*u(inp)
      !V(INP) = perturb/100.
      !W(INP) = perturb/100.



!.....USER INITIALIZATION:

!.....Power law profile:
      !U(INP) = Umag*(zc(inp)/blthick)**expp
      !V(INP) = small
      !W(INP) = small
      !TE(INP)=-3.71354692995576e-5*ZC(INP)**7+0.0008508042*ZC(INP)**6-0.0074501233*ZC(INP)**5 &
      !        +0.0287202493*ZC(INP)**4-0.0279210481*ZC(INP)**3-0.1060775951*ZC(INP)**2 &
      !        +0.1756108394*ZC(INP)+0.2104104465
      !ED(INP)=CMU75*TE(INP)**1.5/(0.07*blthick) !/(CMU*TE(INP))

!.....Bolund hill prescribed inlet:
!      U(IJK)=LOG((ZC(IJK)*120-ZC(INBC)*120)/0.0003) ! log profil case 270
!      V(IJK)=0.
!      W(IJK)=0.
!      TE(IJK)=0.928 ! vrednost data za blind comparison
!      ED(IJK)=(0.4)**3/(CAPPA*(ZC(IJK)-ZC(INBC))*120)  !CMU75*TE(IJK)**1.5/(0.4*0.1) 

!      UO(INP) = small
!      VO(INP) = small
!      WO(INP) = small

!.....log profile for Dobric case;
!      inbc = LI(I)+J 
!      if ( (zc(inp)-zc(inbc)) .lt. 500. ) then
!      !if(k.eq.2) print*, i, j, (zc(inp)-zc(inbc))
!        ust=0.345843
!        U(INP) = min(-(ust/cappa)*log((zc(inp)-zc(inbc))/zzero),-small)
!        TE(INP) = (ust)**2/sqrt(cappa)*(1.-((zc(inp)-zc(inbc))/500.))**2
!      else
!        U(INP) = U(INP)
!        TE(INP) = SMALL
!      endif
!      ed(inp) = max(TE(INP)**1.5/10.,small) !<-- epsilon eq.
!      !ed(inp) = ed(inp)/0.02973/TE(INP) !<-- omega eq.
!      !ed(inp)=ed(inp)/(cmu*te(inp)+small)! k-eps --> k-omega

      ! ! Effective viscosity
      ! IF(LTURB) THEN
      !   IF(STDKEPS.or.DURBIN.or.RNG.or.REALIZABLE.or.LowRE_LB) VIS(INP)=VIS(INP)+DENSIT*TE(INP)**2*CMU/(ED(INP)+SMALL)   
      !   IF(Wilcox.or.SST.or.SAS.or.EARSM_WJ.or.EARSM_M) VIS(INP)=VIS(INP)+DENSIT*TE(INP)/(ED(INP)+SMALL)
      ! ENDIF


      enddo

  !-------------------------------------------------------    
  ! Field initialisation over inner cells + boundary faces
  !-------------------------------------------------------

  ! Density
  den(:) = densit
  ! Effective viscosity
  vis(:)=viscos
  ! ! Temperature
   t(:)=tin
  ! ! Temperature variance
  ! vart(:)=vartin
  ! ! Concentration
  ! con(:)=conin

  ! Reynolds stress tensor components
  uu(:) = 0.0_dp
  vv(:) = 0.0_dp
  ww(:) = 0.0_dp
  uv(:) = 0.0_dp
  uw(:) = 0.0_dp
  vw(:) = 0.0_dp

  ! Turbulent heat fluxes
   utt(:) = 0.0_dp
   vtt(:) = 0.0_dp
   wtt(:) = 0.0_dp

  ! Reynolds stress anisotropy
  ! if(earsm_wj.or.earsm_m) bij(:,:)=0.0_dp

  ! Pressure and pressure correction
  p(:) = small
  pp(:) = p(:)



! 7)  Read Restart File And Set Field Values
  if(lread) then
    call readfiles
    pp(:) = p(:)
  end if



! 8)  Initial Gradient Calculation
  dUdxi=0.0_dp
  dVdxi=0.0_dp
  dWdxi=0.0_dp
  dPdxi=0.0_dp
  dTEdxi=0.0_dp
  dEDdxi=0.0_dp


  if (lstsq .or. lstsq_qr .or. lstsq_dm) then
    call create_lsq_gradients_matrix(U,dUdxi)
  endif

  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)


! 9) Calculate distance to the nearest wall.

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
