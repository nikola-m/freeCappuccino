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
! Face fluxes of velocity for inner faces.
!
!***********************************************************************
! 
  implicit none


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
    ue = face_value(ijp, ijn, xf, yf, zf, fxp, u, dUdxi )
    ve = face_value(ijp, ijn, xf, yf, zf, fxp, v, dVdxi )
    we = face_value(ijp, ijn, xf, yf, zf, fxp, w, dWdxi )
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    ue = face_value(ijn, ijp, xf, yf, zf, fxn, u, dUdxi )
    ve = face_value(ijn, ijp, xf, yf, zf, fxn, v, dVdxi )
    we = face_value(ijn, ijp, xf, yf, zf, fxn, w, dWdxi )
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
    ue = face_value(ijp, ijn, xf, yf, zf, fxp, u, dUdxi )
    ve = face_value(ijp, ijn, xf, yf, zf, fxp, v, dVdxi )
    we = face_value(ijp, ijn, xf, yf, zf, fxp, w, dWdxi )
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    ue = face_value(ijn, ijp, xf, yf, zf, fxn, u, dUdxi )
    ve = face_value(ijn, ijp, xf, yf, zf, fxn, v, dVdxi )
    we = face_value(ijn, ijp, xf, yf, zf, fxn, w, dWdxi )
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
! Facefluxuvw routine used for inlet and outlet boundaries. 
! Constant gradient is assumed between cell center and boundary face center.
!
!***********************************************************************
! 
  implicit none


  integer, intent(in) :: ijp, ijb
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flomass
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
  real(dp) :: fdue,fdve,fdwe,fdui,fdvi,fdwi
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
!     We calculate explicit and implicit diffusion fde and fdi,
!     later we put their difference (fde-fdi) to rhs vector
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Explicit diffussion 
  fdue = game*( (duxii+duxii)*arx + (duyii+dvxii)*ary + (duzii+dwxii)*arz )
  fdve = game*( (duyii+dvxii)*arx + (dvyii+dvyii)*ary + (dvzii+dwyii)*arz )
  fdwe = game*( (duzii+dwxii)*arx + (dwyii+dvzii)*ary + (dwzii+dwzii)*arz )

  ! Implicit diffussion 
  fdui = game*are/dpn*(duxi*xpn+duyi*ypn+duzi*zpn)
  fdvi = game*are/dpn*(dvxi*xpn+dvyi*ypn+dvzi*zpn)
  fdwi = game*are/dpn*(dwxi*xpn+dwyi*ypn+dwzi*zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  ! > Explicit convection: 
  ! - None -
  

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
