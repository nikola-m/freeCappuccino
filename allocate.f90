!***********************************************************************
!
subroutine allocate_arrays
!
!***********************************************************************
!
  use parameters
  use geometry
  use indexes
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


  allocate(x(numNodes),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: x" 
  allocate(y(numNodes),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: y" 
  allocate(z(numNodes),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: z" 

  allocate(xc(numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xc" 
  allocate(yc(numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: yc" 
  allocate(zc(numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: zc" 

  allocate(vol(numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: vol"

  allocate(wallDistance(numCells),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: wallDistance" 


  allocate(arx(numInnerFaces),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: arx" 
  allocate(ary(numInnerFaces),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: ary" 
  allocate(arz(numInnerFaces),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: arz" 

  allocate( xf(numInnerFaces),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xf"
  allocate( yf(numInnerFaces),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: yf"
  allocate( zf(numInnerFaces),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: zf"

  allocate(facint(numInnerFaces),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: facint" 


  allocate(owner(numInnerFaces),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: li"
  allocate(neighbour(numInnerFaces),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: lk"


  allocate( ijl(noc),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: ijl"
  allocate( ijr(noc),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: ijr"



  allocate( xni(ninl),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xni"
  allocate( yni(ninl),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: yni"
  allocate( zni(ninl),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: zni"
  allocate( xfi(ninl),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xfi"
  allocate( yfi(ninl),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: yfi"
  allocate( zfi(ninl),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: zfi"


  allocate( xno(nout),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xno"
  allocate( yno(nout),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: yno"
  allocate( zno(nout),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: zno"
  allocate( xfo(nout),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xfo"
  allocate( yfo(nout),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: yfo"
  allocate( zfo(nout),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: zfo"

  allocate( srds(nsym),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: srds"
  allocate( dns(nsym),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dns"
  allocate( xns(nsym),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xns"
  allocate( yns(nsym),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: yns"
  allocate( zns(nsym),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: zns"
  allocate( xfs(nsym),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xfs"
  allocate( yfs(nsym),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: yfs"
  allocate( zfs(nsym),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: zfs"


  allocate( srdw(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: srdw"
  allocate( dnw(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: dnw"
  allocate( visw(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: visw"
  allocate( ypl(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: ypl"
  allocate( xnw(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xnw"
  allocate( ynw(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: ynw"
  allocate( znw(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: znw"
  allocate( xfw(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xfw"
  allocate( yfw(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: yfw"
  allocate( zfw(nwal),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: zfw"

  allocate( xnpr(npru),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xnpr"
  allocate( ynpr(npru),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: ynpr"
  allocate( znpr(npru),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: znpr"
  allocate( xfpr(npru),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xfpr"
  allocate( yfpr(npru),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: yfpr"
  allocate( zfpr(npru),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: zfpr"

  allocate( srdoc(noc),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: srdoc"
  allocate( xnoc(noc),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xnoc"
  allocate( ynoc(noc),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: ynoc"
  allocate( znoc(noc),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: znoc"
  allocate( xfoc(noc),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: xfoc"
  allocate( yfoc(noc),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: yfoc"
  allocate( zfoc(noc),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: zfoc"
  allocate( foc(noc),stat=ierr) 
    if(ierr /= 0)write(*,*)"allocation error: foc"
 

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
  allocate(den(numCells),stat=ierr) 
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


subroutine deallocate_arrays
!
!     DEALLOCATE ARRAYS
!
  use parameters
  use geometry
  use indexes
  use variables
  use hcoef
  use title_mod
  use statistics
  use omega_turb_models
!
  implicit none

  ! GEOMETRY RELATED

  if (allocated(x)) deallocate(x)
  if (allocated(y)) deallocate(y)
  if (allocated(z)) deallocate(z)

  if (allocated(xc)) deallocate(xc)
  if (allocated(yc)) deallocate(yc)
  if (allocated(zc)) deallocate(zc)

  if (allocated(vol)) deallocate(vol)

  if (allocated(arx)) deallocate(arx)
  if (allocated(ary)) deallocate(ary)
  if (allocated(arz)) deallocate(arz)

  if ( allocated( xf )) deallocate(xf)
  if ( allocated( yf )) deallocate(yf)
  if ( allocated( zf )) deallocate(zf)

  if (allocated(facint)) deallocate(facint)

  if (allocated(wallDistance)) deallocate(wallDistance)

  if (allocated(owner)) deallocate(owner)
  if (allocated(neighbour)) deallocate(neighbour) 

  if ( allocated( ijl )) deallocate(ijl)
  if ( allocated( ijr )) deallocate(ijr)

  if ( allocated( xni )) deallocate(xni)
  if ( allocated( yni )) deallocate(yni)
  if ( allocated( zni )) deallocate(zni)
  if ( allocated( xfi )) deallocate(xfi)
  if ( allocated( yfi )) deallocate(yfi)
  if ( allocated( zfi )) deallocate(zfi)

  if ( allocated( xno )) deallocate(xno)
  if ( allocated( yno )) deallocate(yno)
  if ( allocated( zno )) deallocate(zno)
  if ( allocated( xfo )) deallocate(xfo)
  if ( allocated( yfo )) deallocate(yfo)
  if ( allocated( zfo )) deallocate(zfo)

  if ( allocated( srds )) deallocate(srds)
  if ( allocated( dns )) deallocate(dns)
  if ( allocated( xns )) deallocate(xns)
  if ( allocated( yns )) deallocate(yns)
  if ( allocated( zns )) deallocate(zns)
  if ( allocated( xfs )) deallocate(xfs)
  if ( allocated( yfs )) deallocate(yfs)
  if ( allocated( zfs )) deallocate(zfs)

  if ( allocated( srdw )) deallocate(srdw)
  if ( allocated( dnw )) deallocate(dnw)
  if ( allocated( visw )) deallocate(visw)
  if ( allocated( ypl )) deallocate(ypl)
  if ( allocated( xnw )) deallocate(xnw)
  if ( allocated( ynw )) deallocate(ynw)
  if ( allocated( znw )) deallocate(znw)
  if ( allocated( xfw )) deallocate(xfw)
  if ( allocated( yfw )) deallocate(yfw)
  if ( allocated( zfw )) deallocate(zfw)

  if ( allocated( xnpr)) deallocate(xnpr)
  if ( allocated( ynpr)) deallocate(ynpr)
  if ( allocated( znpr)) deallocate(znpr)
  if ( allocated( xfpr)) deallocate(xfpr)
  if ( allocated( yfpr)) deallocate(yfpr)
  if ( allocated( zfpr)) deallocate(zfpr)

  if ( allocated( srdoc )) deallocate(srdoc)
  if ( allocated( xnoc )) deallocate(xnoc)
  if ( allocated( ynoc )) deallocate(ynoc)
  if ( allocated( znoc )) deallocate(znoc)
  if ( allocated( xfoc )) deallocate(xfoc)
  if ( allocated( yfoc )) deallocate(yfoc)
  if ( allocated( zfoc )) deallocate(zfoc)
  if ( allocated( foc )) deallocate(foc)


  ! VARIABLES
  if (allocated(U)) deallocate(U)
  if (allocated(V)) deallocate(V)
  if (allocated(W)) deallocate(W)
  if (allocated(Flmass)) deallocate(Flmass)
  if (allocated(P)) deallocate(P)
  if (allocated(PP)) deallocate(PP)
  ! if (allocated(TE)) deallocate(TE)
  ! if (allocated(ED)) deallocate(ED)
  ! if (allocated(T)) deallocate(T)
  if (allocated(VIS)) deallocate(VIS)
  if (allocated(DEN)) deallocate(DEN)
  if (allocated(GEN)) deallocate(GEN)
  if (allocated(magStrain)) deallocate(magStrain)
  if (allocated(Vorticity)) deallocate(Vorticity)

  ! if (allocated(VART)) deallocate(VART)
  ! if (allocated(EDD)) deallocate(EDD)
  ! if (allocated(CON)) deallocate(CON)

  ! if (allocated(UTT)) deallocate(UTT)
  ! if (allocated(VTT)) deallocate(VTT)
  ! if (allocated(WTT)) deallocate(WTT)

  ! if (allocated(RET)) deallocate(RET)


  if (allocated(UU)) deallocate(UU)
  if (allocated(UV)) deallocate(UV)
  if (allocated(UW)) deallocate(UW)
  if (allocated(VV)) deallocate(VV)
  if (allocated(VW)) deallocate(VW)
  if (allocated(WW)) deallocate(WW)

  if (allocated(fmi)) deallocate(fmi)
  if (allocated(fmo)) deallocate(fmo)

  ! if (allocated(ALPH)) deallocate(ALPH)
  ! if (allocated(AL_LES)) deallocate(AL_LES)
  ! if (allocated(AL_RANS)) deallocate(AL_RANS)
  ! if (allocated(DIFF)) deallocate(DIFF)

  ! if (allocated(TIMELIMIT)) deallocate(TIMELIMIT)

                
  ! if (allocated(DOMEGA)) deallocate(DOMEGA) 
  ! if (allocated(ALPHASST)) deallocate(ALPHASST)
  ! if (allocated(BETTASST)) deallocate(BETTASST)
  ! if (allocated(PRTINV_TE)) deallocate(PRTINV_TE) 
  ! if (allocated(PRTINV_ED)) deallocate(PRTINV_ED) 

  ! if (allocated(QSAS)) deallocate(QSAS)
  ! if (allocated(LVK)) deallocate(LVK) 
                   
  ! if (allocated(CMUEFF)) deallocate(CMUEFF) 
  ! if (allocated(BIJ)) deallocate(BIJ) 


  ! COEFFICIENT ARRAYS FOR PISO
  if (allocated(h)) deallocate(h)


  ! TIME
  if (allocated(UO)) deallocate(UO)
  if (allocated(VO)) deallocate(VO)
  if (allocated(WO)) deallocate(WO)
  ! if (allocated(TO)) deallocate(TO)
  ! if (allocated(TEO)) deallocate(TEO)
  ! if (allocated(EDO)) deallocate(EDO)
  ! if (allocated(VARTO)) deallocate(VARTO)
  ! if (allocated(CONO)) deallocate(CONO)

  if (allocated(UOO)) deallocate(UOO)
  if (allocated(VOO)) deallocate(VOO)
  if (allocated(WOO)) deallocate(WOO)
  ! if (allocated(TOO)) deallocate(TOO)
  ! if (allocated(TEOO)) deallocate(TEOO)
  ! if (allocated(EDOO)) deallocate(EDOO)
  ! if (allocated(VARTOO)) deallocate(VARTOO)
  ! if (allocated(CONOO)) deallocate(CONOO)


  ! STATISTICS
  if (allocated(U_AVER)) deallocate(U_AVER)
  if (allocated(V_AVER)) deallocate(V_AVER)
  if (allocated(W_AVER)) deallocate(W_AVER)

  if (allocated(TE_AVER)) deallocate(TE_AVER)

  ! if (allocated(T_AVER)) deallocate(T_AVER)

  if (allocated(UU_AVER)) deallocate(UU_AVER)
  if (allocated(VV_AVER)) deallocate(VV_AVER)
  if (allocated(WW_AVER)) deallocate(WW_AVER)
  if (allocated(UV_AVER)) deallocate(UV_AVER)
  if (allocated(UW_AVER)) deallocate(UW_AVER)
  if (allocated(VW_AVER)) deallocate(VW_AVER)
  ! if (allocated(UT_AVER)) deallocate(UT_AVER)
  ! if (allocated(VT_AVER)) deallocate(VT_AVER)
  ! if (allocated(WT_AVER)) deallocate(WT_AVER)
  ! if (allocated(TT_AVER)) deallocate(TT_AVER)


  ! GRADIENT
  if (allocated(dUdxi)) deallocate(dUdxi)
  if (allocated(dVdxi)) deallocate(dVdxi)
  if (allocated(dWdxi)) deallocate(dWdxi)
  ! if (allocated(dTEdxi)) deallocate(dTEdxi)
  ! if (allocated(dEDdxi)) deallocate(dEDdxi)
  if (allocated(dPdxi)) deallocate(dPdxi) 
  ! if (allocated(dTdxi)) deallocate(dTdxi) 
  ! if (allocated(dVartdxi)) deallocate(dVartdxi) 
  ! if (allocated(dCondxi)) deallocate(dCondxi) 
  ! if (allocated(D)) deallocate(D)

end subroutine
