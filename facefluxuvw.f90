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
  use interpolation 

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
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz
  real(dp) :: ixi1,ixi2,ixi3
  real(dp) :: dpn,costheta,costn
  real(dp) :: xi,yi,zi
  real(dp) :: cp,ce
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
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

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
  ! costn = costheta
  ! Orthogonal correction: nrelax =  0 : 
  costn = 1.0_dp
  ! Over-relaxed approach: nrelax = -1 :
  ! costn = 1./costheta
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
  game = vis(ijp)*fxp+vis(ijn)*fxn

  ! Difusion coefficient
  de = game*are/dpn


  ! > Equation coefficients - implicit diffusion and convection
  ce = min(flomass,zero) 
  cp = max(flomass,zero)

  can = -de + min(flomass,zero)
  cap = -de - max(flomass,zero)



  ! > Explicit diffusion: 

  ! Coordinates of point j'
  xi = xc(ijp)*fxp+xc(ijn)*fxn
  yi = yc(ijp)*fxp+yc(ijn)*fxn
  zi = zc(ijp)*fxp+zc(ijn)*fxn


  ! Interpolate gradients defined at cv centers to faces
  duxi = dUdxi(1,ijp)*fxp+dUdxi(1,ijn)*fxn
  duyi = dUdxi(2,ijp)*fxp+dUdxi(2,ijn)*fxn
  duzi = dUdxi(3,ijp)*fxp+dUdxi(3,ijn)*fxn

  ! du/dx_i interpolated at cell face:
  duxii = duxi*d1x + arx/vole*( u(ijn)-u(ijp)-duxi*d2x-duyi*d2y-duzi*d2z ) 
  duyii = duyi*d1y + ary/vole*( u(ijn)-u(ijp)-duxi*d2x-duyi*d2y-duzi*d2z ) 
  duzii = duzi*d1z + arz/vole*( u(ijn)-u(ijp)-duxi*d2x-duyi*d2y-duzi*d2z ) 

  dvxi = dVdxi(1,ijp)*fxp+dVdxi(1,ijn)*fxn
  dvyi = dVdxi(2,ijp)*fxp+dVdxi(2,ijn)*fxn
  dvzi = dVdxi(3,ijp)*fxp+dVdxi(3,ijn)*fxn

  ! dv/dx_i interpolated at cell face:
  dvxii = dvxi*d1x + arx/vole*( v(ijn)-v(ijp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 
  dvyii = dvyi*d1y + ary/vole*( v(ijn)-v(ijp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 
  dvzii = dvzi*d1z + arz/vole*( v(ijn)-v(ijp)-dvxi*d2x-dvyi*d2y-dvzi*d2z ) 

  dwxi = dWdxi(1,ijp)*fxp+dWdxi(1,ijn)*fxn
  dwyi = dWdxi(2,ijp)*fxp+dWdxi(2,ijn)*fxn
  dwzi = dWdxi(3,ijp)*fxp+dWdxi(3,ijn)*fxn

  ! dw/dx_i interpolated at cell face:
  dwxii = dwxi*d1x + arx/vole*( w(ijn)-w(ijp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 
  dwyii = dwyi*d1y + ary/vole*( w(ijn)-w(ijp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 
  dwzii = dwzi*d1z + arz/vole*( w(ijn)-w(ijp)-dwxi*d2x-dwyi*d2y-dwzi*d2z ) 

!---------------------------------------------------------------------------------------
!     We calculate explicit and implicit diffsion fde and fdi,
!     later se put their difference (fde-fdi) to rhs vector:
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

  ! Initialize explicit convective fluxes for higher-order schemes
  fuhigh = 0.0_dp
  fvhigh = 0.0_dp
  fwhigh = 0.0_dp


  ! Explicit convective fluxes for CDS
  if(lcds) then

    ! > Velocities at cell face center

    !  |________uj'_________|_______________ucorr___________________|
    ue=u(ijp)*fxp+u(ijn)*fxn+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))
    !  |________vj'_________|_______________vcorr___________________|
    ve=v(ijp)*fxp+v(ijn)*fxn+(dvxi*(xf-xi)+dvyi*(yf-yi)+dvzi*(zf-zi))
    !  |________wj'_________|_______________wcorr___________________|
    we=w(ijp)*fxp+w(ijn)*fxn+(dwxi*(xf-xi)+dwyi*(yf-yi)+dwzi*(zf-zi))
 

    fuhigh=flomass*ue
    fvhigh=flomass*ve
    fwhigh=flomass*we

  end if


  !
  ! > Convective schemes done as in FLUENT.
  !
  
  if(lcds_flnt) then

    ue = face_value_central(ijp, ijn, xf, yf, zf, u, dUdxi)
    ve = face_value_central(ijp, ijn, xf, yf, zf, v, dVdxi)
    we = face_value_central(ijp, ijn, xf, yf, zf, w, dWdxi)

    fuhigh=flomass*ue
    fvhigh=flomass*ve
    fwhigh=flomass*we

  endif

  if(l2nd_flnt) then

    fuhigh = cp*face_value_2nd_upwind(ijp, xf, yf, zf, u, dUdxi)+&
             ce*face_value_2nd_upwind(ijn, xf, yf, zf, u, dUdxi)
    fvhigh = cp*face_value_2nd_upwind(ijp, xf, yf, zf, v, dVdxi)+&
             ce*face_value_2nd_upwind(ijn, xf, yf, zf, v, dVdxi)
    fwhigh = cp*face_value_2nd_upwind(ijp, xf, yf, zf, w, dWdxi)+&
             ce*face_value_2nd_upwind(ijn, xf, yf, zf, w, dWdxi)

  end if

  if(l2ndlim_flnt) then

    fuhigh = cp*face_value_2nd_upwind_slope_limited(ijp, xf, yf, zf, u, dUdxi, umin, umax)+&
             ce*face_value_2nd_upwind_slope_limited(ijn, xf, yf, zf, u, dUdxi, umin, umax)
    fvhigh = cp*face_value_2nd_upwind_slope_limited(ijp, xf, yf, zf, v, dVdxi, vmin, vmax)+&
             ce*face_value_2nd_upwind_slope_limited(ijn, xf, yf, zf, v, dVdxi, vmin, vmax)
    fwhigh = cp*face_value_2nd_upwind_slope_limited(ijp, xf, yf, zf, w, dWdxi, wmin, wmax)+&
             ce*face_value_2nd_upwind_slope_limited(ijn, xf, yf, zf, w, dWdxi, wmin, wmax)

  end if

  if(lmuscl_flnt) then

    fuhigh = cp*face_value_muscl(ijp, ijn, xf, yf, zf, u, dUdxi)+&
             ce*face_value_muscl(ijn, ijp, xf, yf, zf, u, dUdxi)
    fvhigh = cp*face_value_muscl(ijp, ijn, xf, yf, zf, v, dVdxi)+&
             ce*face_value_muscl(ijn, ijp, xf, yf, zf, v, dVdxi)
    fwhigh = cp*face_value_muscl(ijp, ijn, xf, yf, zf, w, dWdxi)+&
             ce*face_value_muscl(ijn, ijp, xf, yf, zf, w, dWdxi)

  end if


!--------------------------------------------------------------------------------------------
!     BOUNDED HIGH-ORDER CONVECTIVE SCHEMES (Waterson & Deconinck JCP 224 (2007) pp. 182-207)
!--------------------------------------------------------------------------------------------
  if(lluds.or.lsmart.or.lavl.or.lmuscl.or.lumist.or.lgamma) then

  !.... find r's. this is universal for all schemes.
  !.....if flow goes from p to e
  r1 = (2*dUdxi(1,ijp)*xpn + 2*dUdxi(2,ijp)*ypn + 2*dUdxi(3,ijp)*zpn)/(u(ijn)-u(ijp)) - 1.0_dp  
  r2 = (2*dVdxi(1,ijp)*xpn + 2*dVdxi(2,ijp)*ypn + 2*dVdxi(3,ijp)*zpn)/(v(ijn)-v(ijp)) - 1.0_dp 
  r3 = (2*dWdxi(1,ijp)*xpn + 2*dWdxi(2,ijp)*ypn + 2*dWdxi(3,ijp)*zpn)/(w(ijn)-w(ijp)) - 1.0_dp 
  !.....if flow goes from e to p
  r4 = (2*dUdxi(1,ijn)*xpn + 2*dUdxi(2,ijn)*ypn + 2*dUdxi(3,ijn)*zpn)/(u(ijp)-u(ijn)) - 1.0_dp 
  r5 = (2*dVdxi(1,ijn)*xpn + 2*dVdxi(2,ijn)*ypn + 2*dVdxi(3,ijn)*zpn)/(v(ijp)-v(ijn)) - 1.0_dp 
  r6 = (2*dWdxi(1,ijn)*xpn + 2*dWdxi(2,ijn)*ypn + 2*dWdxi(3,ijn)*zpn)/(w(ijp)-w(ijn)) - 1.0_dp  


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

  !=====gamma scheme================================
  elseif(lgamma) then
  !.....psi for gamma scheme:
  !.....if flow goes from p to e
  ! psiw1 = max(0., min(r1, 2.*r1/(r1+1.)))
  ! psiw2 = max(0., min(r2, 2.*r2/(r2+1.)))
  ! psiw3 = max(0., min(r3, 2.*r3/(r3+1.)))
  ! !.....if flow goes from e to p
  ! psie1 = max(0., min(r4, 2.*r4/(r4+1.)))
  ! psie2 = max(0., min(r5, 2.*r5/(r5+1.)))
  ! psie3 = max(0., min(r6, 2.*r6/(r6+1.)))
  !.....psi for koren scheme:
  !.....if flow goes from p to e
       psiw1 = max(0., min(2.*r1, twothirds*r1+onethird, 2.))
       psiw2 = max(0., min(2.*r2, twothirds*r2+onethird, 2.))
       psiw3 = max(0., min(2.*r3, twothirds*r3+onethird, 2.))
  !.....if flow goes from e to p
       psie1 = max(0., min(2.*r4, twothirds*r4+onethird, 2.))
       psie2 = max(0., min(2.*r5, twothirds*r5+onethird, 2.))
       psie3 = max(0., min(2.*r6, twothirds*r6+onethird, 2.))
  !=====end gamma scheme=============================

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
  fuhigh = ce*(u(ijn) + fxn*psie1*(u(ijp)-u(ijn)))+ &
           cp*(u(ijp) + fxp*psiw1*(u(ijn)-u(ijp)))
  !       mass flux| bounded interpolation of velocity to face |

  fvhigh = ce*(v(ijn) + fxn*psie2*(v(ijp)-v(ijn)))+ &
           cp*(v(ijp) + fxp*psiw2*(v(ijn)-v(ijp)))
   
  fwhigh = ce*(w(ijn) + fxn*psie3*(w(ijp)-w(ijn)))+ &
           cp*(w(ijp) + fxp*psiw3*(w(ijn)-w(ijp)))

!.....END OF BOUNDED HIGH-ORDER SCHEMES
!--------------------------------------------------------------------------------------------
  END IF 


! > Explicit part of diffusion fluxes and sources due to deffered correction.

  sup = -gam*(fuhigh-fuuds)+fdue-fdui
  svp = -gam*(fvhigh-fvuds)+fdve-fdvi
  swp = -gam*(fwhigh-fwuds)+fdwe-fdwi

end subroutine




! Some other higher order schemes:
 
  !.....psi for koren scheme:
  !.....if flow goes from p to e
  !      psiw1 = max(0., min(2.*r1, twothirds*r1+onethird, 2.))
  !      psiw2 = max(0., min(2.*r2, twothirds*r2+onethird, 2.))
  !      psiw3 = max(0., min(2.*r3, twothirds*r3+onethird, 2.))
  !.....if flow goes from e to p
  !      psie1 = max(0., min(2.*r4, twothirds*r4+onethird, 2.))
  !      psie2 = max(0., min(2.*r5, twothirds*r5+onethird, 2.))
  !      psie3 = max(0., min(2.*r6, twothirds*r6+onethird, 2.))

  !.....psi for gpl-1/3-alpha-3/2 scheme:  !!!!new scheme>>>koren i ova schema su jako slicne
  !.....if flow goes from p to e
  !      psiw1 = max(0., min(1.5*r1, twothirds*r1+onethird, 2.))
  !      psiw2 = max(0., min(1.5*r2, twothirds*r2+onethird, 2.))
  !      psiw3 = max(0., min(1.5*r3, twothirds*r3+onethird, 2.))
  !.....if flow goes from e to p
  !      psie1 = max(0., min(1.5*r4, twothirds*r4+onethird, 2.))
  !      psie2 = max(0., min(1.5*r5, twothirds*r5+onethird, 2.))
  !      psie3 = max(0., min(1.5*r6, twothirds*r6+onethird, 2.))

  !.....psi for smarter; charm notable; isnas
  !.....if flow goes from p to e
  !      psiw1 = (r1+abs(r1))*(3*r1+1.)/(2*(r1+1.)**2)
  !      psiw2 = (r2+abs(r2))*(3*r2+1.)/(2*(r2+1.)**2)
  !      psiw3 = (r3+abs(r3))*(3*r3+1.)/(2*(r3+1.)**2)
  !.....if flow goes from e to p
  !      psie1 = (r4+abs(r4))*(3*r4+1.)/(2*(r4+1.)**2)
  !      psie2 = (r5+abs(r5))*(3*r5+1.)/(2*(r5+1.)**2)
  !      psie3 = (r6+abs(r6))*(3*r6+1.)/(2*(r6+1.)**2)

  !.....psi for ospre
  !.....if flow goes from p to e
  !      psiw1 = 1.5*r1*(r1+1.)/(r1**2+r1+1.)
  !      psiw2 = 1.5*r2*(r2+1.)/(r2**2+r2+1.)
  !      psiw3 = 1.5*r3*(r3+1.)/(r3**2+r3+1.)
  !.....if flow goes from e to p
  !      psie1 = 1.5*r4*(r4+1.)/(r4**2+r4+1.)
  !      psie2 = 1.5*r5*(r5+1.)/(r5**2+r5+1.)
  !      psie3 = 1.5*r6*(r6+1.)/(r6**2+r6+1.)

  !.....psi for bsou-blui-chakravarthy-osher scheme:
  !.....if flow goes from p to e
  !      psiw1 = max(0., min(2.*r1,1.))
  !      psiw2 = max(0., min(2.*r2,1.))
  !      psiw3 = max(0., min(2.*r3,1.))
  !.....if flow goes from e to p
  !      psie1 = max(0., min(2.*r4,1.))
  !      psie2 = max(0., min(2.*r5,1.))
  !     psie3 = max(0., min(2.*r6,1.))

  !=====spl-3/5 scheme===============================
  !     if(lumist.eq.1) then
  !.....psi for spl-3/5 scheme (new scheme derived from waterson&deconinck's symmetric piecewise-linear scheme):
  !.....if flow goes from p to e
  !     psiw1 = max(0., min(2.*r1, 0.8*r1+0.2, 0.2*r1+0.8, 2.))
  !      psiw2 = max(0., min(2.*r2, 0.8*r2+0.2, 0.2*r2+0.8, 2.))
  !      psiw3 = max(0., min(2.*r3, 0.8*r3+0.2, 0.2*r3+0.8, 2.))
  !.....if flow goes from e to p
  !      psie1 = max(0., min(2.*r4, 0.8*r4+0.2, 0.2*r4+0.8, 2.))
  !      psie2 = max(0., min(2.*r5, 0.8*r5+0.2, 0.2*r5+0.8, 2.))
  !      psie3 = max(0., min(2.*r6, 0.8*r6+0.2, 0.2*r6+0.8, 2.))
  !=====end umist scheme=============================
  !      endif
