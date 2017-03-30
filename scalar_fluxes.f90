module scalar_fluxes
!
! Implementation of common functions for obtaining discretized fluxes for
! transport equations of salars fields. 
! To the outside world we show only the interface function 'facefluxsc', 
! wisely checking the arguments, module decides what function to call.
!
  use types
  use geometry, only: numTotal, numCells

  implicit none


  interface facefluxsc
    module procedure facefluxsc_interior
    module procedure facefluxsc_interior_nonconst_prtr
    module procedure facefluxsc_boundary
  end interface


  private 

  public :: facefluxsc


contains


!***********************************************************************
!
subroutine facefluxsc_interior(ijp, ijn, xf, yf, zf, arx, ary, arz, &
                               flmass, lambda, gam, FI, dFidxi, &
                               prtr, cap, can, suadd, fimin, fimax)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: xc,yc,zc
  use variables, only: vis
  use interpolation

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flmass
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam 
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numCells), intent(in) :: dFidxi
  real(dp), intent(in) :: prtr
  real(dp), intent(inout) :: cap, can, suadd, fimin, fimax


! Local variables
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz,ixi1,ixi2,ixi3,dpn,costheta,costn
  ! real(dp) :: xi,yi,zi
  real(dp) :: Cp,Ce
  real(dp) :: fii,fm
  real(dp) :: fdfie,fdfii,fcfie,fcfii,ffic
  real(dp) :: d1x,d1y,d1z
  real(dp) :: de, vole, game, viste
  real(dp) :: fxp,fxn
  real(dp) :: xpp,ypp,zpp,xep,yep,zep,xpnp,ypnp,zpnp,volep
  real(dp) :: nablaFIxdnnp,nablaFIxdppp
  real(dp) :: dfixi,dfiyi,dfizi
  real(dp) :: dfixii,dfiyii,dfizii
  real(dp) :: r1,r2
  real(dp) :: psie,psiw
!----------------------------------------------------------------------

  dfixi = 0.0_dp
  dfiyi = 0.0_dp
  dfizi = 0.0_dp

  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)


  ! Cell face diffussion coefficient
  viste = (vis(ijp)-viscos)*fxp+(vis(ijn)-viscos)*fxn
  game = (viste*prtr+viscos)


  ! Difusion coefficient for linear system
  de = game*are/dpn

  ! Convection fluxes - uds
  fm = flmass
  ce = min(fm,zero) 
  cp = max(fm,zero)

  !-------------------------------------------------------
  ! System matrix coefficients
  !-------------------------------------------------------
  cap = -de - max(fm,zero)
  can = -de + min(fm,zero)
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Explicit higher order convection
  !-------------------------------------------------------
  if(lcds) then
    !---------------------------------------------
    ! CENTRAL DIFFERENCING SCHEME (CDS) 
    !---------------------------------------------
    ! Interpolate variable FI defined at CV centers to face using corrected CDS:
    ! ! Coordinates of interpolation point j'
    ! xi=xc(ijp)*fxp+xc(ijn)*fxn
    ! yi=yc(ijp)*fxp+yc(ijn)*fxn
    ! zi=zc(ijp)*fxp+zc(ijn)*fxn
    !   |________Ue'___________|_______________Ucorr_____________________|
    fii=fi(ijp)*fxp+fi(ijn)*fxn!+dfixi*(xf-xi)+dfiyi*(yf-yi)+dfizi*(zf-zi)

    ! Explicit second order convection 
    fcfie=fm*fii

  elseif(lluds.or.l2ndlim_flnt.or.l2nd_flnt) then
    !---------------------------------------------
    ! 2ND ORDER UPWIND DIFFERENCING SCHEME (LUDS) 
    !---------------------------------------------
    fcfie = cp*face_value_2nd_upwind_slope_limited(ijp, xf, yf, zf, fi, dFidxi, fimin, fimax)+&
            ce*face_value_2nd_upwind_slope_limited(ijn, xf, yf, zf, fi, dFidxi, fimin, fimax)
  else
    !---------------------------------------------
    ! MUSCL SCHEME (MUSCL)
    !---------------------------------------------
    ! Flux limiter formulation for 'r' coefficients from:
    ! Darwish-Moukalled TVD schemes for unstructured grids, IJHMT, 2003. 
    !---------------------------------------------
    ! Find r's - the gradient ratio. This is universal for all schemes.
    ! If flow goes from P to E
    r1 = (2*dFidxi(1,ijp)*xpn + 2*dFidxi(2,ijp)*ypn + 2*dFidxi(3,ijp)*zpn)/(FI(ijn)-FI(ijp)) - 1.0_dp  
    ! If flow goes from E to P
    r2 = (2*dFidxi(1,ijn)*xpn + 2*dFidxi(2,ijn)*ypn + 2*dFidxi(3,ijn)*zpn)/(FI(ijp)-FI(ijn)) - 1.0_dp 
    ! Find Psi for [ MUSCL ] :
    psiw = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
    psie = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
    ! High order flux at cell face
    fcfie =  ce*(fi(ijn) + fxn*psie*(fi(ijp)-fi(ijn)))+ &
             cp*(fi(ijp) + fxp*psiw*(fi(ijn)-fi(ijp)))
  endif

  !-------------------------------------------------------
  ! Explicit first order convection
  !-------------------------------------------------------
  fcfii = ce*fi(ijn)+cp*fi(ijp)

  !-------------------------------------------------------
  ! Deffered correction for convection = gama_blending*(high-low)
  !-------------------------------------------------------
  ffic = gam*(fcfie-fcfii)


  !-------------------------------------------------------
  ! Explicit part of diffusion
  !-------------------------------------------------------

  ! Components of the unit vector i_ksi
  ixi1=xpn/dpn
  ixi2=ypn/dpn
  ixi3=zpn/dpn

  ! Unit vectors of the face normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are

  ! Angle between vectorsa n and i_xi - we need cosine
  costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

  ! Relaxation factor for higher-order cell face gradient
  ! Minimal correction: nrelax = +1 :
  !costn = costheta
  ! Orthogonal correction: nrelax =  0 : 
  costn = 1.0_dp
  ! Over-relaxed approach: nrelax = -1 :
  !costn = 1./costheta
  ! In general, nrelax can be any signed integer from some 
  ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
  !costn = costheta**nrelax

  ! dpp_j * sf
  vole=xpn*arx+ypn*ary+zpn*arz

  ! Find points P' and Pj'
  xpp=xf-(xf-xc(ijp))*nxx
  ypp=yf-(yf-yc(ijp))*nyy 
  zpp=zf-(zf-zc(ijp))*nzz

  xep=xf-(xf-xc(ijn))*nxx 
  yep=yf-(yf-yc(ijn))*nyy 
  zep=zf-(zf-zc(ijn))*nzz     

  xpnp = xep-xpp 
  ypnp = yep-ypp 
  zpnp = zep-zpp

  volep = arx*xpnp+ary*ypnp+arz*zpnp

 ! Overrelaxed correction vector d2, where S=dpn+d2
  d1x = costn
  d1y = costn
  d1z = costn
  
  xpnp = xpnp*costn
  ypnp = ypnp*costn
  zpnp = zpnp*costn

  ! Interpolate gradients defined at CV centers to faces
  dfixi = dFidxi(1,ijp)*fxp+dFidxi(1,ijn)*fxn
  dfiyi = dFidxi(2,ijp)*fxp+dFidxi(2,ijn)*fxn
  dfizi = dFidxi(3,ijp)*fxp+dFidxi(3,ijn)*fxn

  ! The cell face interpolated gradient (d phi / dx_i)_j:
  ! Nonorthogonal corrections:         ___
  ! nablaFIxdnnp =>> dot_product(dFidxi,dNN')
  ! And:                               ___
  ! nablaFIxdnnp =>> dot_product(dFidxi,dPP')
  nablaFIxdnnp = dFidxi(1,ijn)*(xep-xc(ijn))+dFidxi(2,ijn)*(yep-yc(ijn))+dFidxi(3,ijn)*(zep-zc(ijn))
  nablaFIxdppp = dFidxi(1,ijp)*(xpp-xc(ijp))+dFidxi(2,ijp)*(ypp-yc(ijp))+dFidxi(3,ijp)*(zpp-zc(ijp))

  dfixii = dfixi!*d1x + arx/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
  dfiyii = dfiyi!*d1y + ary/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
  dfizii = dfizi!*d1z + arz/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 


  ! Explicit diffusion
  fdfie = game*(dfixii*arx + dfiyii*ary + dfizii*arz)  

  ! Implicit diffussion 
  fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)


  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  suadd = -ffic+fdfie-fdfii 
  !-------------------------------------------------------

end subroutine



!***********************************************************************
!
subroutine facefluxsc_interior_nonconst_prtr(ijp, ijn, xf, yf, zf, arx, ary, arz, &
                               flmass, lambda, gam, FI, dFidxi, &
                               prtr_ijp,prtr_ijn, cap, can, suadd, fimin, fimax)
!
!***********************************************************************
!
! Models such as k-omega SST have turbulent Prandtl-Schmidt numbers (sigmas)
! not constant but as changing troughout the field. In particular in the SST model,
! the blending function 'fsst' blends two values of sigmas at each cell center.
! The two values are blended and their values at owner and neighbour cells is
! passed into this subroutine to be interpolated to cell face in question.
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: xc,yc,zc
  use variables, only: vis
  use interpolation

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flmass
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam 
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numCells), intent(in) :: dFidxi
  real(dp), intent(in) :: prtr_ijp, prtr_ijn
  real(dp), intent(inout) :: cap, can, suadd, fimin, fimax


! Local variables
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz,ixi1,ixi2,ixi3,dpn,costheta,costn
  ! real(dp) :: xi,yi,zi
  real(dp) :: Cp,Ce
  real(dp) :: fii,fm
  real(dp) :: fdfie,fdfii,fcfie,fcfii,ffic
  real(dp) :: d1x,d1y,d1z
  real(dp) :: de, vole, game, viste, prtr
  real(dp) :: fxp,fxn
  real(dp) :: xpp,ypp,zpp,xep,yep,zep,xpnp,ypnp,zpnp,volep
  real(dp) :: nablaFIxdnnp,nablaFIxdppp
  real(dp) :: dfixi,dfiyi,dfizi
  real(dp) :: dfixii,dfiyii,dfizii
  real(dp) :: r1,r2
  real(dp) :: psie,psiw
!----------------------------------------------------------------------

  dfixi = 0.0_dp
  dfiyi = 0.0_dp
  dfizi = 0.0_dp

  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)


  ! Cell face diffussion coefficient
  viste = (vis(ijp)-viscos)*fxp+(vis(ijn)-viscos)*fxn
  prtr = prtr_ijp*fxp+prtr_ijn*fxn
  game = (viste*prtr+viscos)


  ! Difusion coefficient for linear system
  de = game*are/dpn

  ! Convection fluxes - uds
  fm = flmass
  ce = min(fm,zero) 
  cp = max(fm,zero)

  !-------------------------------------------------------
  ! System matrix coefficients
  !-------------------------------------------------------
  cap = -de - max(fm,zero)
  can = -de + min(fm,zero)
  !-------------------------------------------------------

  !-------------------------------------------------------
  ! Explicit higher order convection
  !-------------------------------------------------------
  if(lcds) then
    !---------------------------------------------
    ! CENTRAL DIFFERENCING SCHEME (CDS) 
    !---------------------------------------------
    ! Interpolate variable FI defined at CV centers to face using corrected CDS:
    ! ! Coordinates of interpolation point j'
    ! xi=xc(ijp)*fxp+xc(ijn)*fxn
    ! yi=yc(ijp)*fxp+yc(ijn)*fxn
    ! zi=zc(ijp)*fxp+zc(ijn)*fxn
    !   |________Ue'___________|_______________Ucorr_____________________|
    fii=fi(ijp)*fxp+fi(ijn)*fxn!+dfixi*(xf-xi)+dfiyi*(yf-yi)+dfizi*(zf-zi)

    ! Explicit second order convection 
    fcfie=fm*fii

  elseif(lluds) then
    !---------------------------------------------
    ! 2ND ORDER UPWIND DIFFERENCING SCHEME (LUDS) 
    !---------------------------------------------
    fcfie = cp*face_value_2nd_upwind_slope_limited(ijp, xf, yf, zf, fi, dFidxi, fimin, fimax)+&
            ce*face_value_2nd_upwind_slope_limited(ijn, xf, yf, zf, fi, dFidxi, fimin, fimax)
  else
    !---------------------------------------------
    ! MUSCL SCHEME (MUSCL)
    !---------------------------------------------
    ! Flux limiter formulation for 'r' coefficients from:
    ! Darwish-Moukalled TVD schemes for unstructured grids, IJHMT, 2003. 
    !---------------------------------------------
    ! Find r's - the gradient ratio. This is universal for all schemes.
    ! If flow goes from P to E
    r1 = (2*dFidxi(1,ijp)*xpn + 2*dFidxi(2,ijp)*ypn + 2*dFidxi(3,ijp)*zpn)/(FI(ijn)-FI(ijp)) - 1.0_dp  
    ! If flow goes from E to P
    r2 = (2*dFidxi(1,ijn)*xpn + 2*dFidxi(2,ijn)*ypn + 2*dFidxi(3,ijn)*zpn)/(FI(ijp)-FI(ijn)) - 1.0_dp 
    ! Find Psi for [ MUSCL ] :
    psiw = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
    psie = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
    ! High order flux at cell face
    fcfie =  ce*(fi(ijn) + fxn*psie*(fi(ijp)-fi(ijn)))+ &
             cp*(fi(ijp) + fxp*psiw*(fi(ijn)-fi(ijp)))
  endif

  !-------------------------------------------------------
  ! Explicit first order convection
  !-------------------------------------------------------
  fcfii = ce*fi(ijn)+cp*fi(ijp)

  !-------------------------------------------------------
  ! Deffered correction for convection = gama_blending*(high-low)
  !-------------------------------------------------------
  ffic = gam*(fcfie-fcfii)


  !-------------------------------------------------------
  ! Explicit part of diffusion
  !-------------------------------------------------------

  ! Components of the unit vector i_ksi
  ixi1=xpn/dpn
  ixi2=ypn/dpn
  ixi3=zpn/dpn

  ! Unit vectors of the face normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are

  ! Angle between vectorsa n and i_xi - we need cosine
  costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

  ! Relaxation factor for higher-order cell face gradient
  ! Minimal correction: nrelax = +1 :
  !costn = costheta
  ! Orthogonal correction: nrelax =  0 : 
  costn = 1.0_dp
  ! Over-relaxed approach: nrelax = -1 :
  !costn = 1./costheta
  ! In general, nrelax can be any signed integer from some 
  ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
  !costn = costheta**nrelax

  ! dpp_j * sf
  vole=xpn*arx+ypn*ary+zpn*arz

  ! Find points P' and Pj'
  xpp=xf-(xf-xc(ijp))*nxx
  ypp=yf-(yf-yc(ijp))*nyy 
  zpp=zf-(zf-zc(ijp))*nzz

  xep=xf-(xf-xc(ijn))*nxx 
  yep=yf-(yf-yc(ijn))*nyy 
  zep=zf-(zf-zc(ijn))*nzz     

  xpnp = xep-xpp 
  ypnp = yep-ypp 
  zpnp = zep-zpp

  volep = arx*xpnp+ary*ypnp+arz*zpnp

 ! Overrelaxed correction vector d2, where S=dpn+d2
  d1x = costn
  d1y = costn
  d1z = costn
  
  xpnp = xpnp*costn
  ypnp = ypnp*costn
  zpnp = zpnp*costn

  ! Interpolate gradients defined at CV centers to faces
  dfixi = dFidxi(1,ijp)*fxp+dFidxi(1,ijn)*fxn
  dfiyi = dFidxi(2,ijp)*fxp+dFidxi(2,ijn)*fxn
  dfizi = dFidxi(3,ijp)*fxp+dFidxi(3,ijn)*fxn

  ! The cell face interpolated gradient (d phi / dx_i)_j:
  ! Nonorthogonal corrections:         ___
  ! nablaFIxdnnp =>> dot_product(dFidxi,dNN')
  ! And:                               ___
  ! nablaFIxdnnp =>> dot_product(dFidxi,dPP')
  nablaFIxdnnp = dFidxi(1,ijn)*(xep-xc(ijn))+dFidxi(2,ijn)*(yep-yc(ijn))+dFidxi(3,ijn)*(zep-zc(ijn))
  nablaFIxdppp = dFidxi(1,ijp)*(xpp-xc(ijp))+dFidxi(2,ijp)*(ypp-yc(ijp))+dFidxi(3,ijp)*(zpp-zc(ijp))

  dfixii = dfixi!*d1x + arx/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
  dfiyii = dfiyi!*d1y + ary/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
  dfizii = dfizi!*d1z + arz/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 


  ! Explicit diffusion
  fdfie = game*(dfixii*arx + dfiyii*ary + dfizii*arz)  

  ! Implicit diffussion 
  fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)


  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  suadd = -ffic+fdfie-fdfii 
  !-------------------------------------------------------

end subroutine



!***********************************************************************
!
subroutine facefluxsc_boundary(ijp, ijn, xf, yf, zf, arx, ary, arz, flmass, FI, dFidxi, prtr, cap, can, suadd)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: xc,yc,zc,numTotal
  use variables, only: vis

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flmass 
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numTotal), intent(in) :: dFidxi
  real(dp), intent(in) :: prtr
  real(dp), intent(inout) :: cap, can, suadd


! Local variables
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz,ixi1,ixi2,ixi3,dpn,costheta,costn
  real(dp) :: Cp,Ce
  real(dp) :: fm

  real(dp) :: fdfie,fdfii!,fcfie,fcfii,ffic

  real(dp) :: d1x,d1y,d1z,d2x,d2y,d2z

  real(dp) :: de, vole, game, viste

  real(dp) :: fxp,fxn
  ! real(dp) :: xpp,ypp,zpp,xep,yep,zep,xpnp,ypnp,zpnp,volep
  ! real(dp) :: nablaFIxdnnp,nablaFIxdppp
  real(dp) :: dfixi,dfiyi,dfizi
  real(dp) :: dfixii,dfiyii,dfizii
  ! real(dp) :: r1,r2
  ! real(dp) :: psie,psiw
!----------------------------------------------------------------------

  dfixi = 0.0_dp
  dfiyi = 0.0_dp
  dfizi = 0.0_dp

  ! > Geometry:

  ! Face interpolation factor
  fxn=1.0_dp
  fxp=0.0_dp

  ! Distance vector between cell center and face center
  xpn=xf-xc(ijp)
  ypn=yf-yc(ijp)
  zpn=zf-zc(ijp)

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! Components of the unit vector i_ksi
  ixi1=xpn/dpn
  ixi2=ypn/dpn
  ixi3=zpn/dpn

  ! Cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! Unit vectors of the normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are

  ! Angle between vectors n and i_xi - we need cosine
  costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

  ! Relaxation factor for higher-order cell face gradient
  ! Minimal correction: nrelax = +1 :
  !costn = costheta
  ! Orthogonal correction: nrelax =  0 : 
  costn = 1.0_dp
  ! Over-relaxed approach: nrelax = -1 :
  !costn = 1./costheta
  ! In general, nrelax can be any signed integer from some 
  ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
  !costn = costheta**nrelax

  ! dpp_j * sf
  vole=xpn*arx+ypn*ary+zpn*arz


  ! Turbulent viscosity
  viste = vis(ijn)-viscos

  ! Cell face diffussion coefficient
  game = viste*prtr+viscos





 !  !-- Intersection point offset and skewness correction --

 !  ! Find points P' and Pj'
 !  xpp=xf-(xf-xc(ijp))*nxx; ypp=yf-(yf-yc(ijp))*nyy; zpp=zf-(zf-zc(ijp))*nzz
 !  xep=xf-(xf-xc(ijn))*nxx; yep=yf-(yf-yc(ijn))*nyy; zep=zf-(zf-zc(ijn))*nzz     

 !  xpnp = xep-xpp; ypnp = yep-ypp; zpnp = zep-zpp
 !  volep = arx*xpnp+ary*ypnp+arz*zpnp

 !  ! Overrelaxed correction vector d2, where S=dpn+d2
 !  d1x = costn
 !  d1y = costn
 !  d1z = costn
  
 !  xpnp = xpnp*costn
 !  ypnp = ypnp*costn
 !  zpnp = zpnp*costn

 !  ! Interpolate gradients defined at CV centers to faces
 !  dfixi = dFidxi(1,ijp)*fxp+dFidxi(1,ijn)*fxn
 !  dfiyi = dFidxi(2,ijp)*fxp+dFidxi(2,ijn)*fxn
 !  dfizi = dFidxi(3,ijp)*fxp+dFidxi(3,ijn)*fxn

 !  ! The cell face interpolated gradient (d phi / dx_i)_j:
 !  ! Nonorthogonal corrections:         ___
 !  ! nablaFIxdnnp =>> dot_product(dFidxi,dNN')
 !  ! And:                               ___
 !  ! nablaFIxdnnp =>> dot_product(dFidxi,dPP')
 !  nablaFIxdnnp = dFidxi(1,ijn)*(xep-xc(ijn))+dFidxi(2,ijn)*(yep-yc(ijn))+dFidxi(3,ijn)*(zep-zc(ijn))
 !  nablaFIxdppp = dFidxi(1,ijp)*(xpp-xc(ijp))+dFidxi(2,ijp)*(ypp-yc(ijp))+dFidxi(3,ijp)*(zpp-zc(ijp))

 !  dfixii = dfixi*d1x + arx/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
 !  dfiyii = dfiyi*d1y + ary/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
 !  dfizii = dfizi*d1z + arz/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 

 !  !-- Intersection point offset and skewness correction --



  !-- Skewness correction --

  ! Overrelaxed correction vector d2, where s=dpn+d2
  d1x = costn
  d1y = costn
  d1z = costn

  d2x = xpn*costn
  d2y = ypn*costn
  d2z = zpn*costn

  ! Interpolate gradients defined at CV centers to faces
  ! dfixi = dFidxi(1,ijp)*fxp+dFidxi(1,ijn)*fxn
  ! dfiyi = dFidxi(2,ijp)*fxp+dFidxi(2,ijn)*fxn
  ! dfizi = dFidxi(3,ijp)*fxp+dFidxi(3,ijn)*fxn
  ! It should be dFidxi(:,ijn), because fxn=1.0, but we write dFidxi(:,ijp) because constant gradient
  ! is applied between cell center and boundary cell face.
  dfixi = dFidxi(1,ijp)
  dfiyi = dFidxi(2,ijp)
  dfizi = dFidxi(3,ijp) 

  !.....du/dx_i interpolated at cell face:
  dfixii = dfixi*d1x + arx/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
  dfiyii = dfiyi*d1y + ary/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
  dfizii = dfizi*d1z + arz/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 

  !-- Skewness correction --
 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Explicit diffusion
  fdfie = game*(dfixii*arx + dfiyii*ary + dfizii*arz)   
  ! Implicit diffussion 
  fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  ! Difusion coefficient
  de = game*are/dpn

  ! Convection fluxes - uds
  fm = flmass
  ce = min(fm,zero) 
  cp = max(fm,zero)

  ! System matrix coefficients
  cap = -de - max(fm,zero)
  can = -de + min(fm,zero)

  ! if(lcds) then
  !   !---------------------------------------------
  !   ! CENTRAL DIFFERENCING SCHEME (CDS) 
  !   !---------------------------------------------
  !   ! Interpolate variable FI defined at CV centers to face using corrected CDS:
  !   !   |________Ue'___________|_______________Ucorr_____________________|
  !   fii=fi(ijp)*fxp+fi(ijn)*fxn!+dfixi*(xf-xi)+dfiyi*(yf-yi)+dfizi*(zf-zi)

  !   ! Explicit second order convection 
  !   fcfie=fm*fii
  ! else
  !   !---------------------------------------------
  !   ! Darwish-Moukalled TVD schemes for unstructured grids, IJHMT, 2003. 
  !   !---------------------------------------------
  !   ! Find r's - the gradient ratio. This is universal for all schemes.
  !   ! If flow goes from P to E
  !   r1 = (2*dFidxi(1,ijp)*xpn + 2*dFidxi(2,ijp)*ypn + 2*dFidxi(3,ijp)*zpn)/(FI(ijn)-FI(ijp)) - 1.0_dp  
  !   ! If flow goes from E to P
  !   r2 = (2*dFidxi(1,ijn)*xpn + 2*dFidxi(2,ijn)*ypn + 2*dFidxi(3,ijn)*zpn)/(FI(ijp)-FI(ijn)) - 1.0_dp 
  !   ! Find Psi for [ MUSCL ] :
  !   psiw = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
  !   psie = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
  !   ! High order flux at cell face
  !   fcfie =  ce*(fi(ijn) + fxn*psie*(fi(ijp)-fi(ijn)))+ &
  !            cp*(fi(ijp) + fxp*psiw*(fi(ijn)-fi(ijp)))
  ! endif

  ! ! Explicit first order convection
  ! fcfii = ce*fi(ijn)+cp*fi(ijp)

  ! Deffered correction for convection = gama_blending*(high-low)
  !ffic = gam*(fcfie-fcfii)

  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  ! suadd = -ffic+fdfie-fdfii 
  suadd = fdfie-fdfii 
  !-------------------------------------------------------

end subroutine


end module