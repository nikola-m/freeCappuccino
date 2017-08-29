module faceflux_velocity
!
! Implementation of common functions for obtaining discretized fluxes for
! Navier-Stokes equation.
! To the outside world we show only the interface function 'facefluxuvw', 
! wisely checking the arguments, module decides what function to call.
!
  use types
  use parameters
  use geometry, only: xc,yc,zc
  use variables
  use gradients, only: sngrad
  use interpolation, only: face_value

  implicit none


  interface facefluxuvw
    module procedure facefluxuvw
    module procedure facefluxuvw_cyclic
    module procedure facefluxuvw_boundary
  end interface


  private 

  public :: facefluxuvw


contains


!***********************************************************************
!
subroutine facefluxuvw(ijp, ijn, xf, yf, zf, arx, ary, arz, flomass, lambda, gam, cap, can, sup, svp, swp)
!
!***********************************************************************
!
  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flomass
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: sup, svp, swp

  ! Local variables
  integer  :: nrelax
  character(len=12) :: approach
  real(dp) :: are,dpn
  real(dp) :: xpn,ypn,zpn
  real(dp) :: cp,ce
  real(dp) :: duxi,duyi,duzi, &
              dvxi,dvyi,dvzi, &
              dwxi,dwyi,dwzi
  real(dp) :: duxii,dvxii,dwxii, &
              duyii,dvyii,dwyii, &
              duzii,dvzii,dwzii
  real(dp) :: de, game
  real(dp) :: fxp,fxn
  real(dp) :: fuuds,fvuds,fwuds,fuhigh,fvhigh,fwhigh
  real(dp) :: ue, ve, we
  real(dp) :: fdue,fdve,fdwe,fdui,fdvi,fdwi
!----------------------------------------------------------------------


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



  ! > Equation coefficients:

  ! Cell face viscosity
  game = vis(ijp)*fxp+vis(ijn)*fxn

  ! Difusion coefficient
  de = game*are/dpn


  ! > Equation coefficients - implicit diffusion and convection
  ce = min(flomass,zero) 
  cp = max(flomass,zero)

  can = -de + min(flomass,zero)
  cap = -de - max(flomass,zero)



  ! > Explicit diffusion: 

  nrelax = 0
  approach  = 'skewness'

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              u, dudxi, nrelax, approach, duxi, duyi, duzi,  &
              duxii, duyii, duzii)

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              v, dvdxi, nrelax, approach, dvxi, dvyi, dvzi, &
              dvxii, dvyii, dvzii)

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              w, dwdxi, nrelax, approach, dwxi, dwyi, dwzi, &
              dwxii, dwyii, dwzii)

!---------------------------------------------------------------------------------------
!     We calculate explicit and implicit diffusion fde and fdi,
!     later we put their difference (fde-fdi) to rhs vector:
!     su = su + (fdue-fdui)
!     sv = sv + (fdve-fdvi)
!     sw = sw + (fdwe-fdwi)
!---------------------------------------------------------------------------------------

  ! Explicit diffussion: 
  fdue = game*( (duxii+duxii)*arx + (duyii+dvxii)*ary + (duzii+dwxii)*arz )
  fdve = game*( (duyii+dvxii)*arx + (dvyii+dvyii)*ary + (dvzii+dwyii)*arz )
  fdwe = game*( (duzii+dwxii)*arx + (dwyii+dvzii)*ary + (dwzii+dwzii)*arz )

  ! Implicit diffussion:
  fdui = game*are/dpn*(duxi*xpn+duyi*ypn+duzi*zpn)
  fdvi = game*are/dpn*(dvxi*xpn+dvyi*ypn+dvzi*zpn)
  fdwi = game*are/dpn*(dwxi*xpn+dwyi*ypn+dwzi*zpn)



  ! > Explicit convection: 

  ! Explicit convective fluxes for UDS
  fuuds = max(flomass,zero)*u(ijp)+min(flomass,zero)*u(ijn)
  fvuds = max(flomass,zero)*v(ijp)+min(flomass,zero)*v(ijn)
  fwuds = max(flomass,zero)*w(ijp)+min(flomass,zero)*w(ijn)


! EXPLICIT CONVECTIVE FLUXES FOR HIGH ORDER BOUNDED SCHEMES
! Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_f[kg/s] * Phi_f[m/s]$

  if( flomass .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    ue = face_value(ijp, ijn, xf, yf, zf, fxp, u, dUdxi, umin, umax)
    ve = face_value(ijp, ijn, xf, yf, zf, fxp, v, dVdxi, vmin, vmax)
    we = face_value(ijp, ijn, xf, yf, zf, fxp, w, dWdxi, wmin, wmax)
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    ue = face_value(ijn, ijp, xf, yf, zf, fxn, u, dUdxi, umin, umax)
    ve = face_value(ijn, ijp, xf, yf, zf, fxn, v, dVdxi, vmin, vmax)
    we = face_value(ijn, ijp, xf, yf, zf, fxn, w, dWdxi, wmin, wmax)
  endif

  fuhigh = flomass*ue
  fvhigh = flomass*ve
  fwhigh = flomass*we



! > Explicit part of diffusion fluxes and sources due to deffered correction.

  sup = -gam*(fuhigh-fuuds)+fdue-fdui
  svp = -gam*(fvhigh-fvuds)+fdve-fdvi
  swp = -gam*(fwhigh-fwuds)+fdwe-fdwi

end subroutine


!***********************************************************************
!
subroutine facefluxuvw_cyclic(ijp, ijn, xf, yf, zf, arx, ary, arz, flomass, lambda, gam, srd, cap, can, sup, svp, swp)
!
!***********************************************************************
!
! Calculation of face fluxes for valocity at domain cut boundaries
! also known as o-c- cuts and also for cyclic domain boundaries.
! What is in common for these two cases is that both cells that share
! the face are actually inner domain cells.
!
!***********************************************************************
!
  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flomass
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam
  real(dp), intent(in) :: srd ! =are/dpn
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: sup, svp, swp

  ! Local variables
  integer  :: nrelax
  character(len=12) :: approach
  real(dp) :: are,dpn
  real(dp) :: xpn,ypn,zpn
  real(dp) :: cp,ce
  real(dp) :: duxi,duyi,duzi, &
              dvxi,dvyi,dvzi, &
              dwxi,dwyi,dwzi
  real(dp) :: duxii,dvxii,dwxii, &
              duyii,dvyii,dwyii, &
              duzii,dvzii,dwzii
  real(dp) :: de, game
  real(dp) :: fxp,fxn
  real(dp) :: fuuds,fvuds,fwuds,fuhigh,fvhigh,fwhigh
  real(dp) :: ue, ve, we
  real(dp) :: fdue,fdve,fdwe,fdui,fdvi,fdwi
!----------------------------------------------------------------------


  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! ! Distance vector between cell centers
  ! xpn=xc(ijn)-xc(ijp)
  ! ypn=yc(ijn)-yc(ijp)
  ! zpn=zc(ijn)-zc(ijp)

  ! ! Distance from P to neighbor N
  ! dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! Distance from P to neighbor N using stored srd (=are/dpn) value
  dpn = are / srd

  xpn = dpn*arx/are
  ypn = dpn*ary/are
  zpn = dpn*arz/are

  ! > Equation coefficients:

  ! Cell face viscosity
  game = vis(ijp)*fxp+vis(ijn)*fxn

  ! Difusion coefficient
  de = game*are/dpn


  ! > Equation coefficients - implicit diffusion and convection
  ce = min(flomass,zero) 
  cp = max(flomass,zero)

  can = -de + min(flomass,zero)
  cap = -de - max(flomass,zero)



  ! > Explicit diffusion: 

  nrelax = 0
  approach  = 'uncorrected'

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              u, dudxi, nrelax, approach, duxi, duyi, duzi,  &
              duxii, duyii, duzii)

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              v, dvdxi, nrelax, approach, dvxi, dvyi, dvzi, &
              dvxii, dvyii, dvzii)

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              w, dwdxi, nrelax, approach, dwxi, dwyi, dwzi, &
              dwxii, dwyii, dwzii)

!---------------------------------------------------------------------------------------
!     We calculate explicit and implicit diffusion fde and fdi,
!     later we put their difference (fde-fdi) to rhs vector:
!     su = su + (fdue-fdui)
!     sv = sv + (fdve-fdvi)
!     sw = sw + (fdwe-fdwi)
!---------------------------------------------------------------------------------------

  ! Explicit diffussion: 
  fdue = game*( (duxii+duxii)*arx + (duyii+dvxii)*ary + (duzii+dwxii)*arz )
  fdve = game*( (duyii+dvxii)*arx + (dvyii+dvyii)*ary + (dvzii+dwyii)*arz )
  fdwe = game*( (duzii+dwxii)*arx + (dwyii+dvzii)*ary + (dwzii+dwzii)*arz )

  ! Implicit diffussion:
  fdui = game*are/dpn*(duxi*xpn+duyi*ypn+duzi*zpn)
  fdvi = game*are/dpn*(dvxi*xpn+dvyi*ypn+dvzi*zpn)
  fdwi = game*are/dpn*(dwxi*xpn+dwyi*ypn+dwzi*zpn)



  ! > Explicit convection: 

  ! Explicit convective fluxes for UDS
  fuuds = max(flomass,zero)*u(ijp)+min(flomass,zero)*u(ijn)
  fvuds = max(flomass,zero)*v(ijp)+min(flomass,zero)*v(ijn)
  fwuds = max(flomass,zero)*w(ijp)+min(flomass,zero)*w(ijn)


! EXPLICIT CONVECTIVE FLUXES FOR HIGH ORDER BOUNDED SCHEMES
! Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_f[kg/s] * Phi_f[m/s]$

  if( flomass .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    ue = face_value(ijp, ijn, xf, yf, zf, fxp, u, dUdxi, umin, umax)
    ve = face_value(ijp, ijn, xf, yf, zf, fxp, v, dVdxi, vmin, vmax)
    we = face_value(ijp, ijn, xf, yf, zf, fxp, w, dWdxi, wmin, wmax)
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    ue = face_value(ijn, ijp, xf, yf, zf, fxn, u, dUdxi, umin, umax)
    ve = face_value(ijn, ijp, xf, yf, zf, fxn, v, dVdxi, vmin, vmax)
    we = face_value(ijn, ijp, xf, yf, zf, fxn, w, dWdxi, wmin, wmax)
  endif

  fuhigh = flomass*ue
  fvhigh = flomass*ve
  fwhigh = flomass*we



! > Explicit part of diffusion fluxes and sources due to deffered correction.

  sup = -gam*(fuhigh-fuuds)+fdue-fdui
  svp = -gam*(fvhigh-fvuds)+fdve-fdvi
  swp = -gam*(fwhigh-fwuds)+fdwe-fdwi

end subroutine



!***********************************************************************
!
subroutine facefluxuvw_boundary(ijp, ijb, xf, yf, zf, arx, ary, arz, flomass, cap, can, sup, svp, swp)
!
!***********************************************************************
!
  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijb
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flomass
  ! real(dp), intent(in) :: lambda
  ! real(dp), intent(in) :: gam
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: sup, svp, swp

! Local variables
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz
  real(dp) :: ixi1,ixi2,ixi3
  real(dp) :: dpn,costheta,costn
  real(dp) :: xi,yi,zi

  real(dp) :: duxi,duyi,duzi, &
                dvxi,dvyi,dvzi, &
                dwxi,dwyi,dwzi

  real(dp) :: duxii,dvxii,dwxii, &
                duyii,dvyii,dwyii, &
                duzii,dvzii,dwzii

  real(dp) :: d2x,d2y,d2z,d1x,d1y,d1z

  real(dp) :: de, vole, game
  real(dp) :: fxp,fxn
  real(dp) :: fuuds,fvuds,fwuds,fuhigh,fvhigh,fwhigh
  real(dp) :: ue, ve, we
  real(dp) :: fdue,fdve,fdwe,fdui,fdvi,fdwi
  real(dp) :: r1,r2,r3,r4,r5,r6
  real(dp) :: psie1,psie2,psie3,psiw1,psiw2,psiw3
!----------------------------------------------------------------------


  ! > Geometry:


  ! Face interpolation factor
  ! fxn=lambda 
  ! fxp=1.0_dp-lambda
  fxn=1.0_dp
  fxp=0.0_dp

  ! Distance vector between cell centers
  ! xpn=xc(ijb)-xc(ijp)
  ! ypn=yc(ijb)-yc(ijp)
  ! zpn=zc(ijb)-zc(ijp)
  xpn=xf-xc(ijp)
  ypn=yf-yc(ijp)
  zpn=zf-zc(ijp)

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! Components of the unit vector i_ksi
  ixi1=xpn/dpn
  ixi2=ypn/dpn
  ixi3=zpn/dpn

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! Unit vectors of the normal
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

  ! dpn * sf
  vole=xpn*arx+ypn*ary+zpn*arz

  ! Overrelaxed correction vector d2, where s=dpn+d2
  d1x = costn
  d1y = costn
  d1z = costn

  d2x = xpn*costn
  d2y = ypn*costn
  d2z = zpn*costn


  ! > Equation coefficients:

  ! Cell face viscosity
  game = vis(ijb) !<- it's (vis(ijp)*fxp+vis(ijb)*fxn) with fxn = 1.0

  ! Difusion coefficient
  de = game*are/dpn


  ! Equation coefficients - implicit diffusion and convection
  can = -de + min(flomass,zero)
  cap = -de - max(flomass,zero)


  ! > Face velocity components and Explicit diffusion: 

  ! Coordinates of point j'
  xi = xf
  yi = yf
  zi = zf


  !.....interpolate gradients defined at cv centers to faces
  ! duxi = dUdxi(1,ijp)*fxp+dUdxi(1,ijb)*fxn
  ! duyi = dUdxi(2,ijp)*fxp+dUdxi(2,ijb)*fxn
  ! duzi = dUdxi(3,ijp)*fxp+dUdxi(3,ijb)*fxn
  duxi = dUdxi(1,ijp)
  duyi = dUdxi(2,ijp)
  duzi = dUdxi(3,ijp) !...because constant gradient

  !.....du/dx_i interpolated at cell face:
  duxii = duxi*d1x + arx/vole*( u(ijb)-u(ijp)-duxi*d2x-duyi*d2y-duzi*d2z ) 
  duyii = duyi*d1y + ary/vole*( u(ijb)-u(ijp)-duxi*d2x-duyi*d2y-duzi*d2z ) 
  duzii = duzi*d1z + arz/vole*( u(ijb)-u(ijp)-duxi*d2x-duyi*d2y-duzi*d2z ) 


  ! dvxi = dVdxi(1,ijp)*fxp+dVdxi(1,ijb)*fxn
  ! dvyi = dVdxi(2,ijp)*fxp+dVdxi(2,ijb)*fxn
  ! dvzi = dVdxi(3,ijp)*fxp+dVdxi(3,ijb)*fxn
  dvxi = dVdxi(1,ijp)
  dvyi = dVdxi(2,ijp)
  dvzi = dVdxi(3,ijp) !...because constant gradient

  !.....dv/dx_i interpolated at cell face:
  dvxii = dvxi*d1x + arx/vole*( v(ijb)-v(ijp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 
  dvyii = dvyi*d1y + ary/vole*( v(ijb)-v(ijp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 
  dvzii = dvzi*d1z + arz/vole*( v(ijb)-v(ijp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 


  ! dwxi = dWdxi(1,ijp)*fxp+dWdxi(1,ijb)*fxn
  ! dwyi = dWdxi(2,ijp)*fxp+dWdxi(2,ijb)*fxn
  ! dwzi = dWdxi(3,ijp)*fxp+dWdxi(3,ijb)*fxn
  dwxi = dWdxi(1,ijp)
  dwyi = dWdxi(2,ijp)
  dwzi = dWdxi(3,ijp) !...because constant gradient

  !.....dw/dx_i interpolated at cell face:
  dwxii = dwxi*d1x + arx/vole*( w(ijb)-w(ijp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 
  dwyii = dwyi*d1y + ary/vole*( w(ijb)-w(ijp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 
  dwzii = dwzi*d1z + arz/vole*( w(ijb)-w(ijp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     we calculate explicit and implicit diffusion fde and fdi,
!     later we put their difference (fde-fdi) to rhs vector
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! explicit diffussion 
  fdue = (duxii+duxii)*arx + (duyii+dvxii)*ary + (duzii+dwxii)*arz
  fdve = (duyii+dvxii)*arx + (dvyii+dvyii)*ary + (dvzii+dwyii)*arz
  fdwe = (duzii+dwxii)*arx + (dwyii+dvzii)*ary + (dwzii+dwzii)*arz

  fdue = game*fdue
  fdve = game*fdve
  fdwe = game*fdwe

  ! implicit diffussion 
  fdui = game*are/dpn*(duxi*xpn+duyi*ypn+duzi*zpn)
  fdvi = game*are/dpn*(dvxi*xpn+dvyi*ypn+dvzi*zpn)
  fdwi = game*are/dpn*(dwxi*xpn+dwyi*ypn+dwzi*zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!++++end: velocities at cell face center and explicit diffusion fluxes+++++++



  ! > Explicit convection: 

  ! Explicit convective fluxes for UDS
  fuuds=max(flomass,zero)*u(ijp)+min(flomass,zero)*u(ijb)
  fvuds=max(flomass,zero)*v(ijp)+min(flomass,zero)*v(ijb)
  fwuds=max(flomass,zero)*w(ijp)+min(flomass,zero)*w(ijb)

  fuhigh=0.0_dp
  fvhigh=0.0_dp
  fwhigh=0.0_dp

  ! Explicit convective fluxes for CDS
  if(lcds) then

    !  |________uj'_________|_______________ucorr___________________|
    ue=u(ijp)*fxp+u(ijb)*fxn+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))

    !  |________vj'_________|_______________vcorr___________________|
    ve=v(ijp)*fxp+v(ijb)*fxn+(dvxi*(xf-xi)+dvyi*(yf-yi)+dvzi*(zf-zi))
    !  |________wj'_________|_______________wcorr___________________|
    we=w(ijp)*fxp+w(ijb)*fxn+(dwxi*(xf-xi)+dwyi*(yf-yi)+dwzi*(zf-zi))
    
    fuhigh=flomass*ue
    fvhigh=flomass*ve
    fwhigh=flomass*we

  end if
!--------------------------------------------------------------------------------------------
!     BOUNDED HIGH-ORDER CONVECTIVE SCHEMES (Waterson & Deconinck JCP 224 (2007) pp. 182-207)
!--------------------------------------------------------------------------------------------
  if(flux_limiter) then

  !.... find r's. this is universal for all schemes.
  !.....if flow goes from p to e
  r1 = (2*dUdxi(1,ijp)*xpn + 2*dUdxi(2,ijp)*ypn + 2*dUdxi(3,ijp)*zpn)/(u(ijb)-u(ijp)) - 1.0_dp  
  r2 = (2*dVdxi(1,ijp)*xpn + 2*dVdxi(2,ijp)*ypn + 2*dVdxi(3,ijp)*zpn)/(v(ijb)-v(ijp)) - 1.0_dp 
  r3 = (2*dWdxi(1,ijp)*xpn + 2*dWdxi(2,ijp)*ypn + 2*dWdxi(3,ijp)*zpn)/(w(ijb)-w(ijp)) - 1.0_dp 
  !.....if flow goes from e to p
  r4 = (2*dUdxi(1,ijb)*xpn + 2*dUdxi(2,ijb)*ypn + 2*dUdxi(3,ijb)*zpn)/(u(ijp)-u(ijb)) - 1.0_dp 
  r5 = (2*dVdxi(1,ijb)*xpn + 2*dVdxi(2,ijb)*ypn + 2*dVdxi(3,ijb)*zpn)/(v(ijp)-v(ijb)) - 1.0_dp 
  r6 = (2*dWdxi(1,ijb)*xpn + 2*dWdxi(2,ijb)*ypn + 2*dWdxi(3,ijb)*zpn)/(w(ijp)-w(ijb)) - 1.0_dp  

  !=====smart scheme================================
  if(lsmart) then
  !.....psi for smart scheme:
  !.....if flow goes from p to e
  psiw1 = max(0., min(2.*r1, 0.75*r1+0.25, 4.))
  psiw2 = max(0., min(2.*r2, 0.75*r2+0.25, 4.))
  psiw3 = max(0., min(2.*r3, 0.75*r3+0.25, 4.))
  !.....if flow goes from e to p
  psie1 = max(0., min(2.*r4, 0.75*r4+0.25, 4.))
  psie2 = max(0., min(2.*r5, 0.75*r5+0.25, 4.))
  psie3 = max(0., min(2.*r6, 0.75*r6+0.25, 4.))
  !=====end smart scheme=============================

  
  !=====avl-smart scheme=============================
  elseif(lavl) then
  !.....psi for avl-smart scheme:
  !.....if flow goes from p to e
  psiw1 = max(0., min(1.5*r1, 0.75*r1+0.25, 2.5))
  psiw2 = max(0., min(1.5*r2, 0.75*r2+0.25, 2.5))
  psiw3 = max(0., min(1.5*r3, 0.75*r3+0.25, 2.5))
  !.....if flow goes from e to p
  psie1 = max(0., min(1.5*r4, 0.75*r4+0.25, 2.5))
  psie2 = max(0., min(1.5*r5, 0.75*r5+0.25, 2.5))
  psie3 = max(0., min(1.5*r6, 0.75*r6+0.25, 2.5))
  !=====end avl-smart scheme==========================
  

  !=====muscl scheme=================================
  elseif(lmuscl) then
  !.....psi for muscl scheme:
  !.....if flow goes from p to e
  psiw1 = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
  psiw2 = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
  psiw3 = max(0., min(2.*r3, 0.5*r3+0.5, 2.))
  !.....if flow goes from e to p
  psie1 = max(0., min(2.*r4, 0.5*r4+0.5, 2.))
  psie2 = max(0., min(2.*r5, 0.5*r5+0.5, 2.))
  psie3 = max(0., min(2.*r6, 0.5*r6+0.5, 2.))
  !=====end muscl scheme=============================
 

  !=====umist scheme=================================
  elseif(lumist) then
  !.....psi for umist scheme:
  !.....if flow goes from p to e
  psiw1 = max(0., min(2.*r1, 0.75*r1+0.25, 0.25*r1+0.75, 2.))
  psiw2 = max(0., min(2.*r2, 0.75*r2+0.25, 0.25*r2+0.75, 2.))
  psiw3 = max(0., min(2.*r3, 0.75*r3+0.25, 0.25*r3+0.75, 2.))
  !.....if flow goes from e to p
  psie1 = max(0., min(2.*r4, 0.75*r4+0.25, 0.25*r4+0.75, 2.))
  psie2 = max(0., min(2.*r5, 0.75*r5+0.25, 0.25*r5+0.75, 2.))
  psie3 = max(0., min(2.*r6, 0.75*r6+0.25, 0.25*r6+0.75, 2.))
  !=====end umist scheme=============================


  elseif(lkoren) then
  ! psi for koren scheme:
  ! if flow goes from p to e
    psiw1 = max(0., min(2.*r1, twothirds*r1+onethird, 2.))
    psiw2 = max(0., min(2.*r2, twothirds*r2+onethird, 2.))
    psiw3 = max(0., min(2.*r3, twothirds*r3+onethird, 2.))
  ! if flow goes from e to p
    psie1 = max(0., min(2.*r4, twothirds*r4+onethird, 2.))
    psie2 = max(0., min(2.*r5, twothirds*r5+onethird, 2.))
    psie3 = max(0., min(2.*r6, twothirds*r6+onethird, 2.))

  elseif(lcharm) then
  ! psi for smarter; charm notable; isnas
  ! if flow goes from p to e
    psiw1 = (r1+abs(r1))*(3*r1+1.)/(2*(r1+1.)**2)
    psiw2 = (r2+abs(r2))*(3*r2+1.)/(2*(r2+1.)**2)
    psiw3 = (r3+abs(r3))*(3*r3+1.)/(2*(r3+1.)**2)
  ! if flow goes from e to p
    psie1 = (r4+abs(r4))*(3*r4+1.)/(2*(r4+1.)**2)
    psie2 = (r5+abs(r5))*(3*r5+1.)/(2*(r5+1.)**2)
    psie3 = (r6+abs(r6))*(3*r6+1.)/(2*(r6+1.)**2)

  elseif(lospre) then
  ! psi for ospre
  ! if flow goes from p to e
    psiw1 = 1.5*r1*(r1+1.)/(r1**2+r1+1.)
    psiw2 = 1.5*r2*(r2+1.)/(r2**2+r2+1.)
    psiw3 = 1.5*r3*(r3+1.)/(r3**2+r3+1.)
  ! if flow goes from e to p
    psie1 = 1.5*r4*(r4+1.)/(r4**2+r4+1.)
    psie2 = 1.5*r5*(r5+1.)/(r5**2+r5+1.)
    psie3 = 1.5*r6*(r6+1.)/(r6**2+r6+1.)

  !=====luds scheme================================
  else
  !.....psi for 2nd order upwind scheme:
  !.....if flow goes from p to e
  psiw1 = 1.0_dp
  psiw2 = 1.0_dp
  psiw3 = 1.0_dp
  !.....if flow goes from e to p
  psie1 = 1.0_dp
  psie2 = 1.0_dp
  psie3 = 1.0_dp
  !=====end luds scheme=============================
  end if



!.....EXPLICIT CONVECTIVE FLUXES FOR HIGH ORDER BOUNDED SCHEMES
!     $Flux_high_order_scheme[N] = mass_flow_rate_in_cell_face_e[kg/s] * Phi_e[m/s]$
!     Phi_e is found by extrapolation from upwind nodes, see eq. (3.29) in Sasa's Thesis.
!     Additional multiplication with PSI is application of flux limiters,
!     see eq. (10) in Waterson&Deconinck paper.


!......Darwish-Moukalled TVD schemes for unstructured girds, IJHMT, 2003.
  fuhigh = min(flomass,zero) * (u(ijb) + fxn*psie1*(u(ijp)-u(ijb)))+ &
           max(flomass,zero) * (u(ijp) + fxp*psiw1*(u(ijb)-u(ijp)))
  !       mass flux| bounded interpolation of velocity to face |

  fvhigh = min(flomass,zero) * (v(ijb) + fxn*psie2*(v(ijp)-v(ijb)))+ &
           max(flomass,zero) * (v(ijp) + fxp*psiw2*(v(ijb)-v(ijp)))
   
  fwhigh = min(flomass,zero) * (w(ijb) + fxn*psie3*(w(ijp)-w(ijb)))+ &
           max(flomass,zero) * (w(ijp) + fxp*psiw3*(w(ijb)-w(ijp)))

!.....END OF BOUNDED HIGH-ORDER SCHEMES
!--------------------------------------------------------------------------------------------
  END IF 


! Explicit part of diffusion fluxes and sources due to deffered correction,
! for all schemes!

  ! sup = -gam*(fuhigh-fuuds)+fdue-fdui
  ! svp = -gam*(fvhigh-fvuds)+fdve-fdvi
  ! swp = -gam*(fwhigh-fwuds)+fdwe-fdwi
  sup = fdue-fdui
  svp = fdve-fdvi
  swp = fdwe-fdwi  !...because gam=0

end subroutine

end module
