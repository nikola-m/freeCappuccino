! ***********************************************************************
!
subroutine init
!
!***********************************************************************
!     Contents:
!
! 1)  Various initialisations
!     1.1)  Parameter Initialisation
!     1.2)  Field Initialisation
! 2)  Read Restart File And Set Field Values
! 3)  Initial Gradient Calculation
! 4)  Initialization of parameters for boundary adjecent cells
! 5)  Calculate distance to the nearest wall.
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use variables
  use title_mod
  use gradients
  use sparse_matrix
  use utils, only: get_unit
  ! use LIS_linear_solver_library
  use output

  implicit none

  ! 
  ! Local variables 
  !
  character(len=80) :: key,field_type,boundary_type
  ! character(len=5)  :: maxno
  ! character(10) :: tol
  integer :: i, ijp, ijn, inp, ini, inw, ijo, ijs, ijb, ioc, iface
  integer :: nfaces,startFace,nFacesOffset
  integer :: input_unit,input_status!,output_unit
  integer :: nsw_backup
  real(dp) :: fxp, fxn, ui, vi, wi
  real(dp) :: nxf, nyf, nzf
  real(dp) :: are
  real(dp) :: sor_backup
  real(dp) :: u0, v0, w0, tke0, ed0

!
!***********************************************************************
!

!
! 1) Various initialisations
!

! 1.1) Parameter Initialisation

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


! 1.2)  Field Initialisation
  
  write(*,'(a)') ' '
  write(*,'(a)') '  Initializing internal field and boundaries (reading 0/.. ):'
  write(*,'(a)') ' '

  ! 
  ! > Velocity


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

          write(*,'(4x,a,3f10.4)') 'uniform',u0,v0,w0

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

      write(*,'(4x,a)') adjustl(boundary_type)

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

  close ( input_unit )

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

  ! O-C- faces
  do i=1,noc
    ioc = iOCStart+i
    ijp = ijl(i)
    ijn = ijr(i)

    fxn = foc(i)
    fxp = 1.0_dp-foc(i)
    
    u(ioc) = u(ijp)*fxp + u(ijn)*fxn
    v(ioc) = v(ijp)*fxp + v(ijn)*fxn
    w(ioc) = w(ijp)*fxp + w(ijn)*fxn
  enddo

  ! ! Pressure Outlet
  ! do i=1,npru
  !   iface = iPressOutletFacesStart + i
  !   ijp = owner(iface)
  !   u(ijp) = ...
  ! enddo
  


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

    close ( input_unit )

    ! What is left are dose boudaries that are not in 0

    ! Symmetry
    do i=1,nsym
      iface = iSymmetryFacesStart + i
      ijp = owner(iface)
      ijs = iSymmetryStart + i
      te(ijs) = te(ijp)
    enddo

    ! OC faces
    do i=1,noc
      ioc = iOCStart+i
      ijp = ijl(i)
      ijn = ijr(i)
      te(ioc) = te(ijp)*(1.0_dp-foc(i)) + te(ijn)*foc(i)
    enddo

    ! ! Pressure Outlet
    ! do i=1,npru
    !   iface = iPressOutletFacesStart + i
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

  close ( input_unit )

  ! What is left are dose boudaries that are not in 0

  ! Symmetry
  do i=1,nsym
    iface = iSymmetryFacesStart + i
    ijp = owner(iface)
    ijs = iSymmetryStart + i
    ed(ijs) = ed(ijp)
  enddo

  ! OC faces
  do i=1,noc
    ioc = iOCStart+i
    ijp = ijl(i)
    ijn = ijr(i)
    ed(ioc) = ed(ijp)*(1.0_dp-foc(i)) + ed(ijn)*foc(i)
  enddo

  ! ! Pressure Outlet
  ! do i=1,npru
  !   iface = iPressOutletFacesStart + i
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

  close ( input_unit )

  ! What is left are dose boudaries that are not in 0

  ! Symmetry
  do i=1,nsym
    iface = iSymmetryFacesStart + i
    ijp = owner(iface)
    ijs = iSymmetryStart + i
    ed(ijs) = ed(ijp)
  enddo

  ! OC faces
  do i=1,noc
    ioc = iOCStart+i
    ijp = ijl(i)
    ijn = ijr(i)
    ed(ioc) = ed(ijp)*(1.0_dp-foc(i)) + ed(ijn)*foc(i)
  enddo

  ! ! Pressure Outlet
  ! do i=1,npru
  !   iface = iPressOutletFacesStart + i
  !   ijp = owner(iface)
  !   ed(ijp) = ...
  ! enddo

  write(*,'(a)') ' '

  endif


  !-------------------------------------------------------    
  ! Field initialisation over inner cells + boundary faces
  !-------------------------------------------------------

  ! Density
  den = densit

  ! Effective viscosity
  vis = viscos

  ! call get_unit ( input_unit )
  ! open ( unit = input_unit, file = '0/vis')
  ! rewind input_unit
  ! do i=1,numCells
  !   read(input_unit,*) vis(i)
  ! enddo
  ! close (input_unit)

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
  if(lturb.and.lasm) bij = 0.0_dp

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

  ! Mass flow for domain (O-C-) and cyclic bounaries
  do i=1,noc
    ioc = iOCStart+i
    ijp = ijl(i)
    ijn = ijr(i)
    iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.

    fmoc(i) = den(ijp)*(arx(iface)*u(ioc)+ary(iface)*v(ioc)+arz(iface)*w(ioc))

  enddo

!
! 2)  Read Restart File And Set Field Values
!
  if(lread) call readfiles


!
! 3)  Initial Gradient Calculation
!
  dUdxi = 0.0_dp
  dVdxi = 0.0_dp
  dWdxi = 0.0_dp
  dPdxi = 0.0_dp
  dTEdxi = 0.0_dp
  dEDdxi = 0.0_dp


! !===================================================
! ! Test gradients:
! u = 0.0_dp
! ! Inlet
! do i=1,ninl
!   iface = iInletFacesStart + i
!   ijp = owner(iface)
!   ijs = iInletStart + i
!   u(ijs) = xf(iface)+yf(iface)+zf(iface)
! enddo
! ! Outlet
! do i=1,nout
!   iface = iOutletFacesStart + i
!   ijp = owner(iface)
!   ijs = iOutletStart + i
!   u(ijs) = xf(iface)+yf(iface)+zf(iface)
! enddo
! ! Symmetry
! do i=1,nsym
!   iface = iSymmetryFacesStart + i
!   ijp = owner(iface)
!   ijs = iSymmetryStart + i
!   u(ijs) = xf(iface)+yf(iface)+zf(iface)
! enddo
! ! Wall
! do i=1,nwal
!   iface = iWallFacesStart + i
!   ijp = owner(iface)
!   ijs = iWallStart + i
!   u(ijs) = xf(iface)+yf(iface)+zf(iface)
! enddo
! u(1:numCells) = xc(1:numCells)+yc(1:numCells)+zc(1:numCells)
! !===================================================


  if (lstsq .or. lstsq_qr .or. lstsq_dm) then
    call create_lsq_gradients_matrix(U,dUdxi)
  endif

  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)
 

! !===================================================
! print*,'gradijenti:'
! do i=1,numCells
!   print*,i,':',dudxi(1,i)
! enddo
! print*,'L0 error norm: ',maxval(abs(1.0d0-dudxi(1,:)))
! print*,'L1 error norm: ',sum(abs(1.0d0-dudxi(1,:)))
! stop
! !===================================================


!
! 4) Calculate distance dnw of wall adjecent cells and distance to the nearest wall of all cell centers.
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


  !
  ! 5) Distance to the nearest wall (needed for some turbulence models).
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
  call laplacian(sv,p) 

  sor_backup = sor(ip)
  nsw_backup = nsw(ip)

  sor(ip) = 1e-10
  nsw(ip) = 500

  ! Solve system
  call iccg(p,ip) 
  ! call bicgstab(p,ip) 
  ! call pmgmres_ilu ( numCells, nnz, ioffset, ja, a, diag, p(1:numCells), ip, su, 100, 4, 1e-8, sor(ip) )
  ! write(maxno,'(i5)') nsw(ip)
  ! write(tol,'(es9.2)') sor(ip)
  ! ! write(options,'(a)') "-i gmres -restart [20] -p ilut -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
  ! write(options,'(a)') "-i cg -p ilu -ilu_fill 1 -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
  ! call solve_csr( numCells, nnz, ioffset, ja, a, su, p )


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


  ! Write wall distance field.
  !+-----------------------------------------------------------------------------+
  ! call get_unit( output_unit )

  ! open(unit=output_unit,file='VTK/wallDistance.vtu')

  ! ! Header
  ! call vtu_write_XML_header ( output_unit )
  ! ! Scalar field
  ! call vtu_write_XML_scalar_field ( output_unit, 'wallDistance', wallDistance )
  ! ! Mesh data
  ! call vtu_write_XML_meshdata ( output_unit )

  ! close( output_unit )
  !+-----------------------------------------------------------------------------+


end subroutine
