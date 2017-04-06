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
! 4)  Various initialisations
!     4.1)  Parameter Initialisation
!     4.2)  Field Initialisation
! 5)  Read Restart File And Set Field Values
! 6)  Initial Gradient Calculation
! 7)  Initialization of parameters for boundary adjecent cells
! 8)  Calculate distance to the nearest wall.
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use variables
  use title_mod
  use gradients
  use sparse_matrix, only: create_CSR_matrix_from_mesh_data,su,sv
  use utils, only: timestamp, show_logo, i4vec_print2, get_unit
  use LIS_linear_solver_library

  implicit none
!
!***********************************************************************
!

! 
! Local variables 
!
  integer :: i, ijp, ijn, inp, ini, inw, ijo, ijs, ijb, iface
  real(dp) :: fxp, fxn, ui, vi, wi
  real(dp) :: nxf, nyf, nzf
  real(dp) :: are

  integer :: nsw_backup
  real(dp) :: sor_backup

  integer :: input_unit,input_status
  integer :: nfaces,startFace,nFacesOffset
  character(80) :: key,field_type,boundary_type
  character(25) :: convective_scheme
  real(dp) :: u0, v0, w0, tke0, ed0

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

!.OPEN & READ INPUT FILE...................................................
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
  READ(5,*) convective_scheme
  READ(5,*) limiter
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
  WRITE(6,'(L1,1x,f4.2,1x,a)') LTRANSIENT,BTIME,'LTRANSIENT,BTIME'
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
  WRITE(6,'(2(L1,1x),1x,a)') BDF,CN,'BDF,CN'
  WRITE(6,'(3(L1,1x),i1,1x,a)') SIMPLE,PISO,PIMPLE,ncorr,'SIMPLE,PISO,PIMPLE,ncorr'
  WRITE(6,'(1(L1,1x),5x,a)') const_mflux,'const_mflux'
  WRITE(6,'(L1,es11.4,5x,a)') CoNumFix, CoNumFixValue,'CoNumFix, CoNumFixValue'
  WRITE(6,'(a)') '---cut here-----------------------------------------------------------------------------'
  WRITE(6,'(a)') ' '


  ! Turbulent flow computation condition
  lturb = levm.or.lasm.or.lles

  ! Switches which define wether we look for some files in folder 0.
  if (lturb) then
    if ( TurbModel==1 .or. TurbModel==2 ) then
      solveEpsilon = .true.
      solveTKE = .true.
    elseif ( TurbModel==3 .or. TurbModel==4 ) then
      solveOmega = .true.
      solveTKE = .true.
    elseif( TurbModel==6 ) then
      solveTKE = .true.
    else
      solveOmega = .false.
      solveEpsilon = .false.
      solveTKE = .false.
    endif
  endif

  ! Set the string for second scalar equation, which appears in linear solver log, default is 'epsilon',
  ! Note, change also the script 'plotResiduals'
  if ( solveOmega ) chvarSolver(6) = 'Omega  '

  ! Choice of convective schemes for velocity

  if(adjustl(convective_scheme) == 'central') then
    lcds = .true.
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
  elseif(adjustl(convective_scheme) == 'gamma') then
    lgamma = .true.
  elseif(adjustl(convective_scheme) == 'central-f') then
    lcds_flnt = .true.
  elseif(adjustl(convective_scheme) == 'linear-f') then
    l2nd_flnt = .true.
  elseif(adjustl(convective_scheme) == 'limited-linear') then
    l2ndlim_flnt = .true.
  elseif(adjustl(convective_scheme) == 'muscl-f') then
    lmuscl_flnt = .true.
  endif

  write(*,'(a)') ' '
  write(*,'(2a)') '  Convective scheme: ', adjustl(convective_scheme)
  write(*,'(a)') ' '

!
! 2)  Open & Read mesh file, calculate mesh geometrical quantities, allocate arrays
!
  call mesh_geometry

  call allocate_arrays

  call allocate_gradients

!
! 3) Index arrays of matrix elements stored in CSR format
!
  call create_CSR_matrix_from_mesh_data


!
! 4) Various initialisations
!


! 4.1) Parameter Initialisation


  ! Reciprocal values of underrelaxation factors
  do i=1,nphi
    urfr(i)=1.0_dp / urf(i)
    urfm(i)=1.0_dp - urf(i)
  enddo

  ! Initial time
  if(.not.lread) time = 0.0_dp

  ! Set to zero cumulative error in continuity
  cumulativeContErr = 0.0_dp

  ! Bulk velocity - important const_mflux flow!
  magUbar = uin


! 4.2)  Field Initialisation
  
  write(*,'(a)') ' '
  write(*,'(a)') '  Initializing internal field and boundaries (reading 0/.. ):'
  write(*,'(a)') ' '

  ! 
  ! > Velocity
  ! 

  call get_unit ( input_unit )
  open ( unit = input_unit, file = '0/U')
  write(*,'(a)') '  0/U'
  rewind input_unit

  do

  read(input_unit,'(a)',iostat = input_status) key

  if(input_status /= 0) exit

    if(adjustl(key)=='internalField') then

      write(*,'(2x,a)') 'internalField'

      read(input_unit,*) field_type

      if(adjustl(field_type)=='uniform') then

          read(input_unit,*) u0,v0,w0

          write(*,'(4x,a,3f9.3)') 'uniform',u0,v0,w0

          do inp = 1,numCells
            u(inp) = u0
            v(inp) = v0
            w(inp) = w0
          enddo

      elseif(adjustl(field_type)=='nonuniform') then

          write(*,'(4x,a)') 'nonuniform'

          do inp = 1,numCells
            read(input_unit,*) u(inp),v(inp),w(inp)
          enddo

      endif

    elseif(adjustl(key)=='boundaryField') then

      write(*,'(2x,a)') 'boundaryField'      
      
      do 

      read(input_unit,'(a)',iostat = input_status) boundary_type

      if(input_status /= 0) exit

      write(*,'(4x,2a)') '>',adjustl(boundary_type)

        if(adjustl(boundary_type) == "inlet") then


            read(input_unit,*) field_type
            if(adjustl(field_type)=='uniform') then

              read(input_unit,*) u0,v0,w0

              write(*,'(6x,a,3f9.3)') 'uniform',u0,v0,w0

              do i = 1,ninl
                iface = iInletFacesStart+i
                ini = iInletStart + i
                u(ini) = u0
                v(ini) = v0
                w(ini) = w0
              enddo     

            else ! 'nonuniform'   

              write(*,'(6x,a)') 'nonuniform'

              do i = 1,ninl
                iface = iInletFacesStart+i
                ini = iInletStart + i
                read(input_unit,*) u(ini),v(ini),w(ini)
              enddo

            endif


        elseif(adjustl(boundary_type) == "outlet") then


            read(input_unit,*) field_type

            if(adjustl(field_type)=='uniform') then

              read(input_unit,*) u0,v0,w0

              write(*,'(6x,a,3f9.3)') 'uniform',u0,v0,w0

              do i=1,nout
                iface = iOutletFacesStart + i
                ijo = iOutletStart+i
                u(ijo) = u0
                v(ijo) = v0
                w(ijo) = w0
              enddo     

            elseif(adjustl(field_type)=='zeroGradient') then  

              write(*,'(6x,a)')  'zeroGradient'

              do i=1,nout
                iface = iOutletFacesStart + i
                ijp = owner(iface)
                ijo = iOutletStart+i
                u(ijo) = u(ijp)
                v(ijo) = v(ijp)
                w(ijo) = w(ijp)
              enddo

            endif



        elseif(adjustl(boundary_type) == "wall") then

            read(input_unit,*) field_type

            if(adjustl(field_type)=='uniform') then ! e.g. moving wall...

              read(input_unit,*) nfaces,startFace
              read(input_unit,*) u0,v0,w0

              write(*,'(6x,a,3f9.3)') 'uniform',u0,v0,w0
              write(*,'(6x,2(a,i8))') 'nfaces:',nfaces,', startFace:',startFace

              nFacesOffset = startFace-iWallFacesStart

              do i = 1,nfaces
                iface = iWallFacesStart + nFacesOffset + i
                inw = iWallStart + nFacesOffset + i
                u(inw) = u0
                v(inw) = v0
                w(inw) = w0
              enddo      

            elseif(adjustl(field_type)=='noSlip') then

              write(*,'(6x,a)') 'noSlip'

              do i = 1,nwal
                iface = iWallFacesStart+i
                inw = iWallStart + i
                u(inw) = zero
                v(inw) = zero
                w(inw) = zero
              enddo

            endif

        endif

      enddo


    endif  

  
  enddo


  ! What is left are dose boudaries that are not in 0

  ! Symmetry
  do i=1,nsym
    iface = iSymmetryFacesStart + i
    ijp = owner(iface)
    ijs = iSymmetryStart + i
    u(ijs) = u(ijp)
    v(ijs) = v(ijp)
    w(ijs) = w(ijp)
  enddo

  ! ! Pressure Outlet
  ! do i=1,npru
  !   iface = iPressOutletFacesStart + i
  !   ijp = owner(iface)
  !   u(ijp) = ...
  ! enddo

  ! ! OC faces
  ! do i=1,noc
  !   iface = iOCFacesStart + i
  !   ijp = owner(iface)
  !   u(ijp) = ...
  ! enddo

  
  ! Initialize minimal and maximal field values for velocity
  umin = minval(u(1:numCells))
  umax = maxval(u(1:numCells))
  vmin = minval(v(1:numCells))
  vmax = maxval(v(1:numCells))
  wmin = minval(w(1:numCells))
  wmax = maxval(w(1:numCells))


  ! 
  ! > TKE Turbulent kinetic energy
  ! 

  if(solveTKE) then

  write(*,'(a)') ' '

  call get_unit ( input_unit )
  open ( unit = input_unit, file = '0/k')
  write(*,'(a)') '  0/k'
  rewind input_unit


  do

  read(input_unit,'(a)',iostat = input_status) key

  if(input_status /= 0) exit

    if(adjustl(key)=='internalField') then

      write(*,'(2x,a)') 'internalField'

      read(input_unit,*) field_type

      if(adjustl(field_type)=='uniform') then

          read(input_unit,*) tke0

          write(*,'(4x,a,f9.3)') 'uniform',tke0

          do inp = 1,numCells
            te(inp) = tke0
          enddo

      elseif(adjustl(field_type)=='nonuniform') then

          write(*,'(4x,a)') 'nonuniform'

          do inp = 1,numCells
            read(input_unit,*) te(inp)
          enddo

      endif

    elseif(adjustl(key)=='boundaryField') then

      write(*,'(2x,a)') 'boundaryField'      
      
      do 

      read(input_unit,'(a)',iostat = input_status) boundary_type

      if(input_status /= 0) exit

      write(*,'(4x,a)') adjustl(boundary_type)

        if(adjustl(boundary_type) == "inlet") then


            read(input_unit,*) field_type
            if(adjustl(field_type)=='uniform') then

              read(input_unit,*) tke0

              write(*,'(6x,a,f9.3)') 'uniform',tke0

              do i = 1,ninl
                iface = iInletFacesStart+i
                ini = iInletStart + i
                te(ini) = tke0
              enddo     

            else ! 'nonuniform'   

              do i = 1,ninl
                iface = iInletFacesStart+i
                ini = iInletStart + i
                read(input_unit,*) te(ini)
              enddo

            endif


        elseif(adjustl(boundary_type) == "outlet") then


            read(input_unit,*) field_type

            if(adjustl(field_type)=='uniform') then

              read(input_unit,*) tke0

              write(*,'(6x,a,f9.3)') 'uniform',tke0

              do i=1,nout
                iface = iOutletFacesStart + i
                ijo = iOutletStart+i
                te(ijo) = tke0
              enddo     

            elseif(adjustl(field_type)=='zeroGradient') then  

              write(*,'(6x,a)')  'zeroGradient'

              do i=1,nout
                iface = iOutletFacesStart + i
                ijp = owner(iface)
                ijo = iOutletStart+i
                te(ijo) = te(ijp)
              enddo

            endif



        elseif(adjustl(boundary_type) == "wall") then

            read(input_unit,*) field_type

            if(adjustl(field_type)=='uniform') then

              read(input_unit,*) tke0

              write(*,'(6x,a,f9.3)') 'uniform',tke0

              do i = 1,nwal
                iface = iWallFacesStart+i
                inw = iWallStart + i
                te(inw) = tke0
              enddo      

            endif

        endif

      enddo


    endif  

  
  enddo

    ! What is left are dose boudaries that are not in 0

    ! Symmetry
    do i=1,nsym
      iface = iSymmetryFacesStart + i
      ijp = owner(iface)
      ijs = iSymmetryStart + i
      te(ijs) = te(ijp)
    enddo

    ! ! Pressure Outlet
    ! do i=1,npru
    !   iface = iPressOutletFacesStart + i
    !   ijp = owner(iface)
    !   te(ijp) = ...
    ! enddo

    ! ! OC faces
    ! do i=1,noc
    !   iface = iOCFacesStart + i
    !   ijp = owner(iface)
    !   te(ijp) = ...
    ! enddo


  endif



  ! 
  ! > ED Turbulent kinetic energy dissipation rate
  ! 

  if(solveEpsilon) then

  write(*,'(a)') ' '

  call get_unit ( input_unit )
  open ( unit = input_unit, file = '0/epsilon')
  write(*,'(a)') '  0/epsilon'
  rewind input_unit


  do

  read(input_unit,'(a)',iostat = input_status) key

  if(input_status /= 0) exit

    if(adjustl(key)=='internalField') then

      write(*,'(2x,a)') 'internalField'

      read(input_unit,*) field_type

      if(adjustl(field_type)=='uniform') then

          read(input_unit,*) ed0

          write(*,'(4x,a,f9.3)') 'uniform',ed0

          do inp = 1,numCells
            ed(inp) = ed0
          enddo

      elseif(adjustl(field_type)=='nonuniform') then

          write(*,'(4x,a)') 'nonuniform'

          do inp = 1,numCells
            read(input_unit,*) ed(inp)
          enddo

      endif

    elseif(adjustl(key)=='boundaryField') then

      write(*,'(2x,a)') 'boundaryField'      
      
      do 

      read(input_unit,'(a)',iostat = input_status) boundary_type

      if(input_status /= 0) exit

      write(*,'(4x,a)') adjustl(boundary_type)

        if(adjustl(boundary_type) == "inlet") then


            read(input_unit,*) field_type
            if(adjustl(field_type)=='uniform') then

              read(input_unit,*) ed0

              write(*,'(6x,a,f9.3)') 'uniform',ed0

              do i = 1,ninl
                iface = iInletFacesStart+i
                ini = iInletStart + i
                ed(ini) = ed0
              enddo     

            else ! 'nonuniform'   

              do i = 1,ninl
                iface = iInletFacesStart+i
                ini = iInletStart + i
                read(input_unit,*) ed(ini)
              enddo

            endif


        elseif(adjustl(boundary_type) == "outlet") then


            read(input_unit,*) field_type

            if(adjustl(field_type)=='uniform') then

              read(input_unit,*) ed0

              write(*,'(6x,a,f9.3)') 'uniform',ed0

              do i=1,nout
                iface = iOutletFacesStart + i
                ijo = iOutletStart+i
                ed(ijo) = ed0
              enddo     

            elseif(adjustl(field_type)=='zeroGradient') then  

              write(*,'(6x,a)')  'zeroGradient'

              do i=1,nout
                iface = iOutletFacesStart + i
                ijp = owner(iface)
                ijo = iOutletStart+i
                ed(ijo) = ed(ijp)
              enddo

            endif



        elseif(adjustl(boundary_type) == "wall") then

            read(input_unit,*) field_type

            if(adjustl(field_type)=='uniform') then

              read(input_unit,*) ed0

              write(*,'(6x,a,f9.3)') 'uniform',ed0

              do i = 1,nwal
                iface = iWallFacesStart+i
                inw = iWallStart + i
                ed(inw) = ed0
              enddo      

            endif

        endif

      enddo


    endif  

  
  enddo

  ! What is left are dose boudaries that are not in 0

  ! Symmetry
  do i=1,nsym
    iface = iSymmetryFacesStart + i
    ijp = owner(iface)
    ijs = iSymmetryStart + i
    ed(ijs) = ed(ijp)
  enddo

  ! ! Pressure Outlet
  ! do i=1,npru
  !   iface = iPressOutletFacesStart + i
  !   ijp = owner(iface)
  !   ed(ijp) = ...
  ! enddo

  ! ! OC faces
  ! do i=1,noc
  !   iface = iOCFacesStart + i
  !   ijp = owner(iface)
  !   ed(ijp) = ...
  ! enddo

  endif



  ! 
  ! > ED Specific turbulent kinetic energy dissipation rate, also turbulence frequency - omega
  ! 

  if(solveOmega) then

  write(*,'(a)') ' '

  call get_unit ( input_unit )
  open ( unit = input_unit, file = '0/omega')
  write(*,'(a)') '  0/omega'
  rewind input_unit

  do

  read(input_unit,'(a)',iostat = input_status) key

  if(input_status /= 0) exit

    if(adjustl(key)=='internalField') then

      write(*,'(2x,a)') 'internalField'

      read(input_unit,*) field_type

      if(adjustl(field_type)=='uniform') then

          read(input_unit,*) ed0

          write(*,'(4x,a,f9.3)') 'uniform',ed0

          do inp = 1,numCells
            ed(inp) = ed0
          enddo

      elseif(adjustl(field_type)=='nonuniform') then

          write(*,'(4x,a)') 'nonuniform'

          do inp = 1,numCells
            read(input_unit,*) ed(inp)
          enddo

      endif

    elseif(adjustl(key)=='boundaryField') then

      write(*,'(2x,a)') 'boundaryField'      
      
      do 

      read(input_unit,'(a)',iostat = input_status) boundary_type

      if(input_status /= 0) exit

      write(*,'(4x,a)') adjustl(boundary_type)

        if(adjustl(boundary_type) == "inlet") then


            read(input_unit,*) field_type
            if(adjustl(field_type)=='uniform') then

              read(input_unit,*) ed0

              write(*,'(6x,a,f9.3)') 'uniform',ed0

              do i = 1,ninl
                iface = iInletFacesStart+i
                ini = iInletStart + i
                ed(ini) = ed0
              enddo     

            else ! 'nonuniform'   

              do i = 1,ninl
                iface = iInletFacesStart+i
                ini = iInletStart + i
                read(input_unit,*) ed(ini)
              enddo

            endif


        elseif(adjustl(boundary_type) == "outlet") then


            read(input_unit,*) field_type

            if(adjustl(field_type)=='uniform') then

              read(input_unit,*) ed0

              write(*,'(6x,a,f9.3)') 'uniform',ed0

              do i=1,nout
                iface = iOutletFacesStart + i
                ijo = iOutletStart+i
                ed(ijo) = ed0
              enddo     

            elseif(adjustl(field_type)=='zeroGradient') then  

              write(*,'(6x,a)')  'zeroGradient'

              do i=1,nout
                iface = iOutletFacesStart + i
                ijp = owner(iface)
                ijo = iOutletStart+i
                ed(ijo) = ed(ijp)
              enddo

            endif



        elseif(adjustl(boundary_type) == "wall") then

            read(input_unit,*) field_type

            if(adjustl(field_type)=='uniform') then

              read(input_unit,*) ed0

              write(*,'(6x,a,f9.3)') 'uniform',ed0

              do i = 1,nwal
                iface = iWallFacesStart+i
                inw = iWallStart + i
                ed(inw) = ed0
              enddo      

            endif

        endif

      enddo


    endif  

  
  enddo

  ! What is left are dose boudaries that are not in 0

  ! Symmetry
  do i=1,nsym
    iface = iSymmetryFacesStart + i
    ijp = owner(iface)
    ijs = iSymmetryStart + i
    ed(ijs) = ed(ijp)
  enddo

  ! ! Pressure Outlet
  ! do i=1,npru
  !   iface = iPressOutletFacesStart + i
  !   ijp = owner(iface)
  !   ed(ijp) = ...
  ! enddo

  ! ! OC faces
  ! do i=1,noc
  !   iface = iOCFacesStart + i
  !   ijp = owner(iface)
  !   ed(ijp) = ...
  ! enddo

  write(*,'(a)') ' '

  endif


! ! Field initialisation loop over inner cells--------------------------------
!   do inp = 1,numCells

! ! Initialization of field variables from input file:
!   u(inp) = uin
!   v(inp) = vin
!   w(inp) = win
!   te(inp) = tein
!   ed(inp) = edin

!   ! Channel flow:
!   ! Random number based fluctuation of mean profile            
!   ! call init_random_seed()
!   ! call random_number(perturb)    
!   ! perturb = 0.9+perturb/5. ! Max perturbation is +/- 10% of mean profile
!   ! u(inp) = perturb*u(inp)
!   ! v(inp) = perturb/100.
!   ! w(inp) = perturb/100.

!   enddo
 
  ! Initialize variables at boundaries

  ! ! Inlet 

  ! ! Outlet - zero gradient
  ! do i=1,nout
  ! iface = iOutletFacesStart + i
  ! ijp = owner(iface)
  ! ijo = iOutletStart+i
  ! u(ijo) = u(ijp)
  ! v(ijo) = v(ijp)
  ! w(ijo) = w(ijp)
  ! te(ijo) = te(ijp)
  ! ed(ijo) = ed(ijp)
  ! enddo

  ! ! Symmetry
  ! do i=1,nsym
  ! iface = iSymmetryFacesStart + i
  ! ijp = owner(iface)
  ! ijs = iSymmetryStart + i
  !     ! Gradient test:
  !     ! u(ijs) = xf(iface)+yf(iface)+zf(iface)
  ! u(ijs) = u(ijp)
  ! v(ijs) = v(ijp)
  ! w(ijs) = w(ijp)
  ! te(ijs) = te(ijp)
  ! ed(ijs) = ed(ijp)
  ! enddo

  ! ! Wall
  ! do i=1,nwal
  ! iface = iWallFacesStart + i
  ! ijw = iWallStart+i
  !     ! Gradient test:
  !     ! u(ijw) = xf(iface)+yf(iface)+zf(iface)
  ! u(ijw) = 0.0_dp
  ! v(ijw) = 0.0_dp
  ! w(ijw) = 0.0_dp
  ! te(ijw) = 0.0_dp
  ! ed(ijw) = 0.0_dp
  ! enddo

  ! Moving Wall
  ! do i=1,20!nwalm
  ! iface = iWallFacesStart + i !iWallMFacesStart + i
  ! inw = iWallStart+i !iWallMStart+i
  ! u(inw) = 1.0_dp
  ! v(inw) = 0.0_dp
  ! w(inw) = 0.0_dp
  ! te(inw) = tein
  ! ed(inw) = edin
  ! enddo

  ! ! Pressure Outlet
  ! do i=1,npru
  ! iface = iPressOutletFacesStart + i
  ! ijp = owner(iface)
  ! u(ijp) = ...
  ! enddo

  ! ! OC faces
  ! do i=1,noc
  ! iface = iOCFacesStart + i
  ! ijp = owner(iface)
  ! u(ijp) = ...
  ! enddo




  !-------------------------------------------------------    
  ! Field initialisation over inner cells + boundary faces
  !-------------------------------------------------------

  ! Density
  den = densit

  ! Effective viscosity
  vis = viscos

  ! Temperature
  if(lcal(ien)) t = tin

  ! Temperature variance
  if(lcal(ivart)) vart = vartin
  
  ! Concentration
  if(lcal(icon)) con=conin

  ! Reynolds stress tensor components
  if (lturb) then
    uu = 0.0_dp
    vv = 0.0_dp
    ww = 0.0_dp
    uv = 0.0_dp
    uw = 0.0_dp
    vw = 0.0_dp
  endif

  ! Turbulent heat fluxes
  if(lcal(ien).and.lbuoy) then
    utt = 0.0_dp
    vtt = 0.0_dp
    wtt = 0.0_dp
  endif

  ! Reynolds stress anisotropy
  if(lasm) bij = 0.0_dp

  ! Pressure and pressure correction
  p = 0.0_dp
  pp = p

  ! Initialize mass flow
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)

    fxn = facint(i)
    fxp = 1.0_dp-facint(i)

    ui = u(ijp)*fxp + u(ijn)*fxn
    vi = v(ijp)*fxp + v(ijn)*fxn
    wi = w(ijp)*fxp + w(ijn)*fxn

    flmass(i) = den(ijp)*(arx(i)*ui+ary(i)*vi+arz(i)*wi)

  enddo


!
! 5)  Read Restart File And Set Field Values
!
  if(lread) then
    call readfiles
    pp = p
  end if


!
! 6)  Initial Gradient Calculation
!
  dUdxi = 0.0_dp
  dVdxi = 0.0_dp
  dWdxi = 0.0_dp
  dPdxi = 0.0_dp
  dTEdxi = 0.0_dp
  dEDdxi = 0.0_dp


  if (lstsq .or. lstsq_qr .or. lstsq_dm) then
    call create_lsq_gradients_matrix(U,dUdxi)
  endif

  ! Gradient limiter:
  write(*,'(a)') ' '

  if(adjustl(limiter) == 'Barth-Jespersen') then
    write(*,*) ' Gradient limiter: Barth-Jespersen'
  elseif(adjustl(limiter) == 'Venkatakrishnan') then
    write(*,*) ' Gradient limiter: Venkatakrishnan'
  elseif(adjustl(limiter) == 'MVenkatakrishnan') then
    write(*,*) ' Gradient limiter: Wang modified Venkatakrishnan'
  else
    write(*,*) ' Gradient limiter: none'
  endif

  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)
 


! print*,'gradijenti:'
! do i=1,numCells
!   print*,i,':',dudxi(1,i)
! enddo
! print*,'L0 error norm: ',maxval(abs(1.0d0-dudxi(1,:)))
! print*,'L1 error norm: ',sum(abs(1.0d0-dudxi(1,:)))
! stop

!
! 7) Calculate distance dnw of wall adjecent cells and distance to the nearest wall of all cell centers.
!

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

    visw(i) = viscos
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
  ! 8) Distance to the nearest wall (needed for some turbulence models).
  !

  write(*,*) ' '
  write(*,*) ' Calculate distance to the nearest wall:'
  write(*,*) ' '

  ! Source term
  su(1:numCells) = -Vol(1:numCells)

  ! Initialize solution
  p = 0.0_dp

  do i=1,numBoundaryFaces
    iface = numInnerFaces+i
    ijp = numCells+i
    p(ijp) = 0.0_dp
  enddo

  ! Wall
  do i=1,nwal
    iface = iWallFacesStart + i
    ijp = iWallStart+i
    p(ijp) = 0.0_dp
  enddo

  !  Coefficient array for Laplacian
  sv = 1.0_dp       

  ! Laplacian operator and BCs         
  call fvm_laplacian(sv,p) 

  sor_backup = sor(ip)
  nsw_backup = nsw(ip)

  sor(ip) = 1e-16
  nsw(ip) = 1000

  ! Solve system
  call iccg(p,ip) 
  ! call gaussSeidel(p,ip) 
  ! call bicgstab(p,ip) 
  ! call dpcg(p,ip)
  ! call solve_csr(numCells,nnz,ioffset,ja,a,su,p)

  ! Update values at constant gradient bc faces - we need these values for correct gradients

  ! Inlet faces
  do i=1,ninl
    iface = iInletFacesStart+i
    ijp = owner(iface)
    ijb = iInletStart+i
    p(ijb)=p(ijp)
  end do

  ! Outlet faces
  do i=1,nout
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i
    p(ijb)=p(ijp)
  end do

  ! Symmetry faces
  do i=1,nsym
    iface = iSymmetryFacesStart+i
    ijp = owner(iface)
    ijb = iSymmetryStart+i
    p(ijb)=p(ijp)
  end do

  sor(ip) = sor_backup
  nsw(ip) = nsw_backup

  ! Gradient of solution field stored in p (gradient stored in dPdxi) :
  call grad(p,dPdxi)

  ! Wall distance computation from Poisson eq. solution stored in pp:
  wallDistance = -sqrt(  dPdxi(1,:)*dPdxi(1,:)+dPdxi(2,:)*dPdxi(2,:)+dPdxi(3,:)*dPdxi(3,:)  ) + &
                  sqrt(  dPdxi(1,:)*dPdxi(1,:)+dPdxi(2,:)*dPdxi(2,:)+dPdxi(3,:)*dPdxi(3,:) + 2*p(1:numCells)  )

  ! Clear arrays
  su = 0.0_dp
  sv = 0.0_dp 
  p = 0.0_dp
  dPdxi = 0.0_dp

  ! write(*,'(a)') ' '
  ! do i=1,numCells
  !   write(6,'(es11.4)') wallDistance(i)
  ! enddo   
  ! stop

end subroutine
