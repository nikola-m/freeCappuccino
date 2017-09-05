!***********************************************************************
!
subroutine facefluxuvw(ijp, ijn, xf, yf, zf, arx, ary, arz, flomass, lambda, gam, cap, can, sup, svp, swp)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: xc,yc,zc
  use variables
  use gradients, only: sngrad
  use interpolation, only: face_value

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
  character(len=8) :: approach
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