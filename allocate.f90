!***********************************************************************
!
subroutine allocate_arrays
!
!***********************************************************************
!
  use parameters
  use geometry, only: numTotal,numCells,numInnerFaces,nnz,ninl,nout,nwal,noc
  use variables
  use hcoef
  use title_mod
  use statistics

  implicit none 
!
!***********************************************************************
!
  integer :: ierr 

  ! Variables

  ! Velocities 
  allocate(u(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: u" 
  allocate(v(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: v" 
  allocate(w(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: w" 

  allocate(uo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: uo" 
  allocate(vo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: vo" 
  allocate(wo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: wo" 

  if( bdf .and. btime.gt.0.99 ) then
    allocate(uoo(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: uoo" 
    allocate(voo(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: voo" 
    allocate(woo(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: woo" 
  endif

  allocate(dUdxi(3,numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dUdxi"
  allocate(dVdxi(3,numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dVdxi"
  allocate(dWdxi(3,numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dWdxi"

  if (ltransient) then
    allocate(u_aver(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: u_aver" 
    allocate(v_aver(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: v_aver" 
    allocate(w_aver(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: w_aver" 
  endif



  ! Pressure and pressure correction
  allocate(p(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: p" 
  allocate(pp(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: pp" 
  allocate(dPdxi(3,numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dPdxi"



  ! Turbulent K.E. and dissipation or any other turbulent scalar taking its place
  allocate(te(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: te" 
  allocate(ed(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: ed" 

  allocate(teo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: teo" 
  allocate(edo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: edo" 

  if( bdf .and. btime.gt.0.99 ) then
    allocate(teoo(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: teoo" 
    allocate(edoo(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: edoo" 
  endif

  allocate(dTEdxi(3,numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dTEdxi"
  allocate(dEDdxi(3,numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dEDdxi"

  if (ltransient) then
    allocate(te_aver(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: te_aver" 
  endif


  ! Effective viscosity  
  allocate(vis(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: vis" 



  ! Temperature
  if(lcal(ien)) then

    allocate(t(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: t"

    allocate(to(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: to"

    if( bdf .and. btime.gt.0.99 ) then
      allocate(too(numTotal),stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: too" 
    endif

    allocate(dTdxi(3,numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: dTdxi"

    if (ltransient) then
      allocate(t_aver(numTotal),stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: t_aver" 
    endif

  endif

  ! Concentration
  if(lcal(icon)) then

    allocate(con(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: con"

    allocate(cono(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: cono"

    if( bdf .and. btime.gt.0.99 ) then
      allocate(conoo(numTotal),stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: conoo" 
    endif

    allocate(dCondxi(3,numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: dCondxi"

  endif

  ! Temperature variance and dissipation of temperature variance
  if(lcal(ivart)) then

    allocate(vart(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: vart" 

    allocate(varto(numTotal),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: varto" 

    if( bdf .and. btime.gt.0.99 ) then
      allocate(vartoo(numTotal),stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: vartoo" 
    endif

    allocate(dVartdxi(3,numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: dVartdxi"

    if (ltransient) then
      allocate(tt_aver(numTotal),stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: tt_aver" 
    endif

    ! allocate(edd(numTotal),stat=ierr) 
    !   if(ierr /= 0)write(*,*)"allocation error: edd"

  endif


  ! Turbulent heat fluxes

  if(lcal(ien).and.lbuoy) then

    allocate(utt(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: utt" 
    allocate(vtt(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: vtt" 
    allocate(wtt(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: wtt"

      if (ltransient) then

        allocate(ut_aver(numCells),stat=ierr) 
          if(ierr /= 0)write(*,*)"allocation error: ut_aver" 
        allocate(vt_aver(numCells),stat=ierr) 
          if(ierr /= 0)write(*,*)"allocation error: vt_aver" 
        allocate(wt_aver(numCells),stat=ierr) 
          if(ierr /= 0)write(*,*)"allocation error: wt_aver"

      endif

  endif

  ! Density
  allocate(den(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: den" 

  ! Mass flows trough east, north and top faces
  allocate(flmass(numInnerFaces),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: flmass" 

  ! Mass flows trough inlet, outlet, and O- C- grid cuts
  allocate(fmi(ninl),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: fmi"
  allocate(fmo(nout),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: fmo"
  allocate(fmoc(noc),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: fmoc"

  allocate( visw(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: visw"
  allocate( ypl(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: ypl"

  ! Turbulence production
  allocate(gen(numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: gen"  

  allocate(magStrain(numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: magStrain" 
  allocate(Vorticity(numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: Vorticity" 



  ! Re_t used in low-re turbulent models
  ! allocate(ret(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ret" 


  ! Reynolds stresses
  if(lturb) then

    allocate(uu(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: uu" 
    allocate(uv(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: uv" 
    allocate(uw(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: uw" 
    allocate(vv(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: vv" 
    allocate(vw(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: vw" 
    allocate(ww(numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: ww" 

    if(lasm) then
      allocate(bij(5,numCells), stat=ierr)
        if(ierr/=0)write(*,*)"allocate error: bij"
     endif

    if(ltransient) then

      allocate(uu_aver(numCells),stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: uu_aver" 
      allocate(vv_aver(numCells),stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: vv_aver" 
      allocate(ww_aver(numCells),stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: ww_aver" 
      allocate(uv_aver(numCells),stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: uv_aver" 
      allocate(uw_aver(numCells),stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: uw_aver" 
      allocate(vw_aver(numCells),stat=ierr) 
        if(ierr /= 0)write(*,*)"allocation error: vw_aver" 

    endif

  endif 


  ! Coefficient arrays for PISO and PIMPLE
  if( piso .or. pimple ) then
    allocate(h(nnz),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: h" 
  endif


end subroutine
