!***********************************************************************
!
subroutine allocate_arrays
!
!***********************************************************************
!
  use parameters
  use geometry
  use variables
  use hcoef
  use title_mod
  use statistics
  use omega_turb_models

  implicit none 
!
!***********************************************************************
!
  integer :: ierr


  ! allocate(x(numNodes),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: x" 
  ! allocate(y(numNodes),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: y" 
  ! allocate(z(numNodes),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: z" 

  ! allocate(xc(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xc" 
  ! allocate(yc(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: yc" 
  ! allocate(zc(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: zc" 

  ! allocate(vol(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: vol"

  ! allocate(wallDistance(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: wallDistance" 


  ! allocate(arx(numInnerFaces),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: arx" 
  ! allocate(ary(numInnerFaces),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ary" 
  ! allocate(arz(numInnerFaces),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: arz" 

  ! allocate( xf(numInnerFaces),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xf"
  ! allocate( yf(numInnerFaces),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: yf"
  ! allocate( zf(numInnerFaces),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: zf"

  ! allocate(facint(numInnerFaces),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: facint" 


  ! allocate(owner(numInnerFaces),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: li"
  ! allocate(neighbour(numInnerFaces),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: lk"


  ! allocate( ijl(noc),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ijl"
  ! allocate( ijr(noc),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ijr"



  ! allocate( xni(ninl),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xni"
  ! allocate( yni(ninl),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: yni"
  ! allocate( zni(ninl),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: zni"
  ! allocate( xfi(ninl),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xfi"
  ! allocate( yfi(ninl),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: yfi"
  ! allocate( zfi(ninl),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: zfi"


  ! allocate( xno(nout),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xno"
  ! allocate( yno(nout),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: yno"
  ! allocate( zno(nout),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: zno"
  ! allocate( xfo(nout),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xfo"
  ! allocate( yfo(nout),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: yfo"
  ! allocate( zfo(nout),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: zfo"

  ! allocate( srds(nsym),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: srds"
  ! allocate( dns(nsym),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: dns"
  ! allocate( xns(nsym),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xns"
  ! allocate( yns(nsym),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: yns"
  ! allocate( zns(nsym),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: zns"
  ! allocate( xfs(nsym),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xfs"
  ! allocate( yfs(nsym),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: yfs"
  ! allocate( zfs(nsym),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: zfs"


  ! allocate( srdw(nwal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: srdw"
  ! allocate( dnw(nwal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: dnw"
  ! allocate( visw(nwal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: visw"
  ! allocate( ypl(nwal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ypl"
  ! allocate( xnw(nwal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xnw"
  ! allocate( ynw(nwal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ynw"
  ! allocate( znw(nwal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: znw"
  ! allocate( xfw(nwal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xfw"
  ! allocate( yfw(nwal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: yfw"
  ! allocate( zfw(nwal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: zfw"

  ! allocate( xnpr(npru),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xnpr"
  ! allocate( ynpr(npru),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ynpr"
  ! allocate( znpr(npru),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: znpr"
  ! allocate( xfpr(npru),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xfpr"
  ! allocate( yfpr(npru),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: yfpr"
  ! allocate( zfpr(npru),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: zfpr"

  ! allocate( srdoc(noc),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: srdoc"
  ! allocate( xnoc(noc),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xnoc"
  ! allocate( ynoc(noc),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ynoc"
  ! allocate( znoc(noc),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: znoc"
  ! allocate( xfoc(noc),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: xfoc"
  ! allocate( yfoc(noc),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: yfoc"
  ! allocate( zfoc(noc),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: zfoc"
  ! allocate( foc(noc),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: foc"
 

  ! Variables

  ! Velocities 
  allocate(u(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: u" 
  allocate(v(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: v" 
  allocate(w(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: w" 

  ! Pressure and pressure correction
  allocate(p(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: p" 
  allocate(pp(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: pp" 

  ! Turbulent K.E. and dissipation or any other turbulent scalar taking its place
  ! allocate(te(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: te" 
  ! allocate(ed(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ed" 

  ! Effective viscosity  
  allocate(vis(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: vis" 

  ! Temperature
  ! allocate(t(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: t"

  ! Concentration
  ! allocate(con(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: con" 

  ! Temperature variance and ..
  ! allocate(vart(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: vart" 
  ! allocate(edd(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: edd"
  ! Density

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


  ! Turbulence production
  allocate(gen(numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: gen"  

  allocate(magStrain(numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: magStrain" 
  allocate(Vorticity(numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: Vorticity" 

  ! Turbulent heat fluxes
  ! allocate(utt(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: utt" 
  ! allocate(vtt(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: vtt" 
  ! allocate(wtt(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: wtt"

  ! Re_t used in low-re turbulent models
  ! allocate(ret(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ret" 


  ! Reynolds stresses
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
 

  ! if (alphamodel) then
  !   allocate(alph(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: alph" 
  !   allocate(al_les(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: al_les" 
  !   allocate(al_rans(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: al_rans" 
  !   allocate(diff(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: diff" 
  ! endif



  ! ! Durbin time-scale limiter
  ! if (durbin) then
  !   allocate(timelimit(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: timelimit" 
  ! endif

  !     for menter sst, sas, earsm_wj (wallin-johansson) and earsm_m (menter implementation of wj)
  ! if (sst.or.sas.or.earsm_wj.or.earsm_m) then 
                    
  !   allocate(domega(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: domega" 
  !   allocate(alphasst(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: alphasst" 
  !   allocate(bettasst(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: bettasst" 
  !   allocate(prtinv_te(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: prtinv_te" 
  !   allocate(prtinv_ed(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: prtinv_ed" 

  ! if (sas) then
  !   allocate(qsas(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: qsas" 
  !   allocate(lvk(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: lvk"
  ! endif

  ! if (earsm_wj.or.earsm_m) then                         
  !   allocate(cmueff(numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: cmueff" 
  !   allocate(bij(5,numCells),stat=ierr) 
  !     if(ierr /= 0)write(*,*)"allocation error: bij"  
  ! endif

  ! endif


  !    Coefficient arrays for PISO and PIMPLE
  if( piso .or. pimple ) then
    allocate(h(nnz),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: h" 
  endif


  !     Time
  allocate(uo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: uo" 
  allocate(vo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: vo" 
  allocate(wo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: wo" 
    
  ! allocate(to(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: to" 
  ! allocate(teo(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: teo" 
  ! allocate(edo(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: edo" 
  ! allocate(varto(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: varto" 
  ! allocate(cono(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: cono" 

  if( bdf .and. btime.gt.0.5 ) then
  allocate(uoo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: uoo" 
  allocate(voo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: voo" 
  allocate(woo(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: woo" 
  ! allocate(too(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: too" 
  ! allocate(teoo(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: teoo" 
  ! allocate(edoo(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: edoo" 
  ! allocate(vartoo(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: vartoo" 
  ! allocate(conoo(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: conoo" 
  endif




  !     Statistics
  if (ltransient) then

  allocate(u_aver(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: u_aver" 
  allocate(v_aver(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: v_aver" 
  allocate(w_aver(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: w_aver" 
  allocate(te_aver(numTotal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: te_aver" 

  ! allocate(t_aver(numTotal),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: t_aver" 

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

  ! allocate(ut_aver(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: ut_aver" 
  ! allocate(vt_aver(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: vt_aver" 
  ! allocate(wt_aver(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: wt_aver" 
  ! allocate(tt_aver(numCells),stat=ierr) 
  !   if(ierr /= 0)write(*,*)"allocation error: tt_aver"

  endif 


  ! Gradient
  allocate(dUdxi(3,numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dUdxi"
  allocate(dVdxi(3,numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dVdxi"
  allocate(dWdxi(3,numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dWdxi"
  allocate(dPdxi(3,numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dPdxi"

      ! allocate(dTEdxi(3,numCells),stat=ierr) 
      !   if(ierr /= 0)write(*,*)"allocation error: dTEdxi"
      ! allocate(dEDdxi(3,numCells),stat=ierr) 
      !   if(ierr /= 0)write(*,*)"allocation error: dEDdxi"


      !allocate(dTdxi(3,numCells),stat=ierr) 
        !if(ierr /= 0)write(*,*)"allocation error: dTdxi"
      !allocate(dVartdxi(3,numCells),stat=ierr) 
        !if(ierr /= 0)write(*,*)"allocation error: dVartdxi"
      !allocate(dCondxi(3,numCells),stat=ierr) 
        !if(ierr /= 0)write(*,*)"allocation error: dCondxi"

      ! !allocate(d(6,numCells),stat=ierr) 
      ! allocate(d(3,6,numCells),stat=ierr) 
      !   if(ierr /= 0)write(*,*)"allocation error: d"

end subroutine
