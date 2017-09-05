subroutine fvm_div(phi,u)
!  
!******************************************************************************
!
!     Fills matrix with coefficients representing implicit FVM discretization 
!     of dicergence operator: div(rho*u).
!
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix

  implicit none

  real(dp), dimension(...), intent(in) :: phi
  real(dp), dimension(numTotal), intent(in) :: u

  !
  ! Local variables
  !

  integer :: i, k, ijp, ijn, ijb, iface
  real(dp) :: cap, can
  real(dp) :: are,dpw
  real(dp) :: gam
  real(dp) :: sup, svp, swp


  ! Initialize matrix array
  a = 0.0_dp

  ! > Assemble Laplacian system matrix

  ! Internal faces:
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxconvection(ijp, ijn, xf(i), yf(i), zf(i), phi(i), facint(i), gam, cap, can, sup, svp, swp)

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_value_index(i)
    a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_value_index(i)
    a(k) = cap

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - can

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - cap

    ! > Sources: 

    su(ijp) = su(ijp) + sup
    sv(ijp) = sv(ijp) + svp
    sw(ijp) = sw(ijp) + swp

    su(ijn) = su(ijn) - sup
    sv(ijn) = sv(ijn) - svp
    sw(ijn) = sw(ijn) - swp

  end do


  ! o- and c-grid cuts
  do i=1,noc

    iface = iOCFacesStart+i
    ijp=ijl(i)
    ijn=ijr(i)

    call facefluxdivergence(ijp, ijn, xf(iface), yf(iface), zf(iface), phi(iface), foc(i), gam, al(i), ar(i), sup, svp, swp)
    
    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - ar(i)

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - al(i)

  end do


!.....Modify matrix coefficients to reflect presence of Boundary Conditions in PDE problem.

  ! Contribution from inlet boundaries
  do i=1,ninl
    iface = iInletFacesStart+i
    ijp = owner(iface)
    ijb = iInletStart+i

    k=diag(ijp)
    are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)
    dpw = sqrt( (xc(ijp)-xf(iface))**2 + (yc(ijp)-yf(iface))**2 + (zc(ijp)-zf(iface))**2 )

    a(k) = a(k) - mu(ijp)*are/dpw !..or mu_wall*are/dpw;  
    su(ijp) = su(ijp) + a(k)*phi(ijb)

  end do

  ! Contribution from outlet boundaries
  do i=1,nout
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i

    k=diag(ijp)
    are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)
    dpw = sqrt( (xc(ijp)-xf(iface))**2 + (yc(ijp)-yf(iface))**2 + (zc(ijp)-zf(iface))**2 )

    a(k) = a(k) - mu(ijp)*are/dpw !..or mu_wall*are/dpw;  
    su(ijp) = su(ijp) + a(k)*phi(ijb)

  end do

  ! Contribution from symmetry boundaries
  do i=1,nsym
    iface = iSymmetryFacesStart
    ijp = owner(iface)
    ijb = iSymmetryStart+i

    k=diag(ijp)

    are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)
    dpw = sqrt( (xc(ijp)-xf(iface))**2 + (yc(ijp)-yf(iface))**2 + (zc(ijp)-zf(iface))**2 )

    a(k) = a(k) - mu(ijp)*are/dpw !..or mu_wall*are/dpw;  

    ! a(k) = a(k) - mu(ijp)*srds(i)  
    su(ijp) = su(ijp) + a(k)*phi(ijb)

  end do

  ! Contribution from wall boundaries
  do i=1,nwal
    iface = iWallFacesStart+i
    ijp = owner(iface)
    ijb = iWallStart+i

    k=diag(ijp)

    are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)
    dpw = sqrt( (xc(ijp)-xf(iface))**2 + (yc(ijp)-yf(iface))**2 + (zc(ijp)-zf(iface))**2 )

    a(k) = a(k) - mu(ijp)*are/dpw !..or mu_wall*are/dpw;  
    ! a(k) = a(k) - mu(ijp)*srdw(i)  
    su(ijp) = su(ijp) + a(k)*phi(ijb)
    
  end do

  ! Contribution from pressure outlet boundaries
  do i=1,npru
    ijp = owner(iPressOutletFacesStart+i)
    ijb = iPressOutletStart+i
  end do


end subroutine





!***********************************************************************
!
subroutine facefluxconvection(ijp, ijn, xf, yf, zf, flomass, lambda, gam, cap, can, sup, svp, swp)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: xc,yc,zc
  use variables

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  ! real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flomass
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: sup, svp, swp

! Local variables
  ! real(dp) :: are
  ! real(dp) :: onethird, twothirds
  real(dp) :: xpn,ypn,zpn
  ! real(dp) :: nxx,nyy,nzz
  ! real(dp) :: ixi1,ixi2,ixi3
  ! real(dp) :: dpn,costheta,costn
  real(dp) :: xi,yi,zi
  real(dp) :: cp,ce

  real(dp) :: duxi,duyi,duzi,dvxi,dvyi,dvzi,dwxi,dwyi,dwzi

  ! real(dp) :: duxii,dvxii,dwxii, &
  !               duyii,dvyii,dwyii, &
  !               duzii,dvzii,dwzii

  ! real(dp) :: d2x,d2y,d2z,d1x,d1y,d1z

  ! real(dp) :: de, vole, game
  real(dp) :: fxp,fxn
  real(dp) :: fuuds,fvuds,fwuds,fuhigh,fvhigh,fwhigh
  real(dp) :: ue, ve, we
  ! real(dp) :: fdue,fdve,fdwe,fdui,fdvi,fdwi
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

  ! Coordinates of point j'
  xi = xc(ijp)*fxp+xc(ijn)*fxn
  yi = yc(ijp)*fxp+yc(ijn)*fxn
  zi = zc(ijp)*fxp+zc(ijn)*fxn

  ! > Equation coefficients:

  can =  min(flomass,zero)
  cap = -max(flomass,zero)


  ! > Explicit part due to deffered correction:


  ! Interpolate gradients defined at cv centers to faces
  duxi = dUdxi(1,ijp)*fxp+dUdxi(1,ijn)*fxn
  duyi = dUdxi(2,ijp)*fxp+dUdxi(2,ijn)*fxn
  duzi = dUdxi(3,ijp)*fxp+dUdxi(3,ijn)*fxn

  dvxi = dVdxi(1,ijp)*fxp+dVdxi(1,ijn)*fxn
  dvyi = dVdxi(2,ijp)*fxp+dVdxi(2,ijn)*fxn
  dvzi = dVdxi(3,ijp)*fxp+dVdxi(3,ijn)*fxn

  dwxi = dWdxi(1,ijp)*fxp+dWdxi(1,ijn)*fxn
  dwyi = dWdxi(2,ijp)*fxp+dWdxi(2,ijn)*fxn
  dwzi = dWdxi(3,ijp)*fxp+dWdxi(3,ijn)*fxn



  ! > Explicit convection: 

  ! Explicit convective fluxes for UDS
  fuuds=max(flomass,zero)*u(ijp)+min(flomass,zero)*u(ijn)
  fvuds=max(flomass,zero)*v(ijp)+min(flomass,zero)*v(ijn)
  fwuds=max(flomass,zero)*w(ijp)+min(flomass,zero)*w(ijn)


  ! Initialize explicit convective fluxes for higher-order schemes
  fuhigh=0.0_dp
  fvhigh=0.0_dp
  fwhigh=0.0_dp

  ! Explicit convective fluxes for CDS
  if(lcds) then

    ! > Velocities at cell face center

    !  |________uj'_________|_______________ucorr___________________|
    ue=u(ijp)*fxp+u(ijn)*fxn+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))
    !  |________vj'_________|_______________vcorr___________________|
    ve=v(ijp)*fxp+v(ijn)*fxn+(dvxi*(xf-xi)+dvyi*(yf-yi)+dvzi*(zf-zi))
    !  |________wj'_________|_______________wcorr___________________|
    we=w(ijp)*fxp+w(ijn)*fxn+(dwxi*(xf-xi)+dwyi*(yf-yi)+dwzi*(zf-zi))

    ! ue = face_interpolated(u,dUdxi,inp,idew,idns,idtb,fxp,fxe)
    ! ve = face_interpolated(v,dVdxi,inp,idew,idns,idtb,fxp,fxe)
    ! we = face_interpolated(w,dWdxi,inp,idew,idns,idtb,fxp,fxe)

    fuhigh=flomass*ue
    fvhigh=flomass*ve
    fwhigh=flomass*we

  end if

!--------------------------------------------------------------------------------------------
!     BOUNDED HIGH-ORDER CONVECTIVE SCHEMES (Waterson & Deconinck JCP 224 (2007) pp. 182-207)
!--------------------------------------------------------------------------------------------
  if(lsmart.or.lavl.or.lmuscl.or.lumist.or.lgamma) then

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
  psiw1 = max(0., min(r1, 2.*r1/(r1+1.)))
  psiw2 = max(0., min(r2, 2.*r2/(r2+1.)))
  psiw3 = max(0., min(r3, 2.*r3/(r3+1.)))
  !.....if flow goes from e to p
  psie1 = max(0., min(r4, 2.*r4/(r4+1.)))
  psie2 = max(0., min(r5, 2.*r5/(r5+1.)))
  psie3 = max(0., min(r6, 2.*r6/(r6+1.)))
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

  ce = min(flomass,zero) 
  cp = max(flomass,zero)

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


! > Explicit part due to deffered correction.

  sup = -gam*(fuhigh-fuuds)
  svp = -gam*(fvhigh-fvuds)
  swp = -gam*(fwhigh-fwuds)


end subroutine



