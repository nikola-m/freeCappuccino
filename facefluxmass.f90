!***********************************************************************
!
subroutine facefluxmass(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, cap, can, flmass)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: xc,yc,zc,vol
  use variables, only: den,U,V,W,dUdxi,dVdxi,dWdxi,p,dpdxi
  use sparse_matrix, only: apu,apv,apw
  use gradients

  implicit none

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: cap, can
  real(dp), intent(inout) :: flmass

  ! Local variables
  real(dp) :: fxn, fxp
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn,dene,smdpn,sfdpnr
  real(dp) :: xi,yi,zi
  real(dp) :: nxx,nyy,nzz,xpp,ypp,zpp,xep,yep,zep
  real(dp) :: ui,vi,wi,ue,ve,we,dpe,dpex,dpey,dpez
  real(dp) :: dpxi,dpyi,dpzi
  real(dp) :: duxi,duyi,duzi


  ! Tentative (!) velocity gradients: 
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)


  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.-lambda

  ! Coordinates of point e'
  xi=xc(ijp)*fxp+xc(ijn)*fxn
  yi=yc(ijp)*fxp+yc(ijn)*fxn
  zi=zc(ijp)*fxp+zc(ijn)*fxn

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! Unit vectors of the normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are


  sfdpnr=1./(arx*xpn*nxx+ary*ypn*nyy+arz*zpn*nzz)

  ! density at the cell face
  dene=den(ijp)*fxp+den(ijn)*fxn

  ! COEFFICIENTS OF PRESSURE-CORRECTION EQUATION
  smdpn = (arx**2+ary**2+arz**2)*sfdpnr
  cap = -(fxp*vol(ijp)*apu(ijp)+fxn*vol(ijn)*apu(ijn))*dene*smdpn
  can = cap

!
! CELL FACE PRESSURE GRADIENTS AND VELOCITIES
!
!////////////////////////////////////////////////////////
!     RHIE-CHOW velcity interolation at face
!
!   Uf=UI+DPDXI-API*Sf*(Pn-Pp)
!         _
!   UI-> (U)f -> second order interpolation at face
!            _______________ 
!   DPDXI-> (dPdx*Vol*(1/ap))f -> second order interpolation at cell face
!          ______
!   API*Sf*(Pn-Pp) -> (1/ap)f*Sf*(Pn-Pp) cell face coefficient 1/Ap x Area_f x (p1-p2)
!  
!   Finally:     
!         __     _______________     ______
!   Uf = (U)f + (dPdx*Vol*(1/ap))f - (1/ap)f*Sf*(Pn-Pp)
!
!   and:
!   Flmass = Densit*dot(Uf,Sf)
!/////////////////////////////////////////////////////////


  ! UI-> (U)f -> second order interpolation at face
  !+Interpolate velocities to face center:+++++++++++++++++++++++++
  ! Interpolate gradients defined at CV centers to faces
  duxi = dudxi(1,ijp)*fxp+dudxi(1,ijn)*fxn
  duyi = dudxi(2,ijp)*fxp+dudxi(2,ijn)*fxn
  duzi = dudxi(3,ijp)*fxp+dudxi(3,ijn)*fxn
  !    |________Ue'_________|_______________Ucorr___________________|
  ui = u(ijp)*fxp+u(ijn)*fxn+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))
  ! UI = face_interpolated(U,dUdxi,ijp,idew,idns,idtb,fxp,fxn)

  duxi = dvdxi(1,ijp)*fxp+dvdxi(1,ijn)*fxn
  duyi = dvdxi(2,ijp)*fxp+dvdxi(2,ijn)*fxn
  duzi = dvdxi(3,ijp)*fxp+dvdxi(3,ijn)*fxn
  !  |________Ve'_________|_______________Vcorr___________________|
  vi=v(ijp)*fxp+v(ijn)*fxn+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))
  ! VI = face_interpolated(V,dVdxi,ijp,idew,idns,idtb,fxp,fxn)

  duxi = dwdxi(1,ijp)*fxp+dwdxi(1,ijn)*fxn
  duyi = dwdxi(2,ijp)*fxp+dwdxi(2,ijn)*fxn
  duzi = dwdxi(3,ijp)*fxp+dwdxi(3,ijn)*fxn
  !  |________We'_________|_______________Wcorr___________________|
  wi=w(ijp)*fxp+w(ijn)*fxn+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)) 
  ! WI = face_interpolated(W,dWdxi,ijp,idew,idns,idtb,fxp,fxn) 
  
  !+END: Interpolate velocities to face center:+++++++++++++++++++++++++


  ! DPDXI-> (dPdx*Vol*(1/ap))f -> second order interpolation at cell face
  !+Interpolate pressure gradients to cell face center++++++++++++++++++
  dpxi = (fxn*Vol(ijn)*Apu(ijn)*dPdxi(1,ijn)+fxp*Vol(ijp)*Apu(ijp)*dPdxi(1,ijp))*xpn*nxx
  dpyi = (fxn*Vol(ijn)*Apv(ijn)*dPdxi(2,ijn)+fxp*Vol(ijp)*Apv(ijp)*dPdxi(2,ijp))*ypn*nyy
  dpzi = (fxn*Vol(ijn)*Apw(ijn)*dPdxi(3,ijn)+fxp*Vol(ijp)*Apw(ijp)*dPdxi(3,ijp))*zpn*nzz
  !+END: Interpolate pressure gradients to cell face center+++++++++++++


  ! (1/ap)f*Sf*(Pn-Pp)
  !+Pressure deriv. along normal+++++++++++++++++++++++++++++++++++++++++ 
  !.....Values at points p' and e' due to non-orthogonality. 
  xpp=xf-(xf-xc(ijp))*nxx; ypp=yf-(yf-yc(ijp))*nyy; zpp=zf-(zf-zc(ijp))*nzz
  xep=xf-(xf-xc(ijn))*nxx; yep=yf-(yf-yc(ijn))*nyy; zep=zf-(zf-zc(ijn))*nzz
  !.....Distances |P'P| and |E'E| projected ionto x,y,z-axis
  xpp=xpp-xc(ijp); ypp=ypp-yc(ijp); zpp=zpp-zc(ijp)
  xep=xep-xc(ijn); yep=yep-yc(ijn); zep=zep-zc(ijn)

  dpe = (p(ijn)-p(ijp)) + &
  ( dPdxi(1,ijn)*xep+dPdxi(2,ijn)*yep+dPdxi(3,ijn)*zep - & !<<--Correction
    dPdxi(1,ijp)*xpp+dPdxi(2,ijp)*ypp+dPdxi(3,ijp)*zpp  )  !<<|

  dpex = (fxn*Vol(ijn)*Apu(ijn)+fxp*Vol(ijp)*Apu(ijp))*(arx*sfdpnr)*dpe
  dpey = (fxn*Vol(ijn)*Apv(ijn)+fxp*Vol(ijp)*Apv(ijp))*(ary*sfdpnr)*dpe
  dpez = (fxn*Vol(ijn)*Apw(ijn)+fxp*Vol(ijp)*Apw(ijp))*(arz*sfdpnr)*dpe
  !+END: Pressure deriv. along normal++++++++++++++++++++++++++++++++++++

  ! Rhie-Chow Interpolation 
  ue = ui - dpex + dpxi
  ve = vi - dpey + dpyi
  we = wi - dpez + dpzi

  ! MASS FLUX via Rhie-Chow Interpolation of velocity
  flmass=dene*(ue*arx+ve*ary+we*arz)

  ! dpxii = (fxn*dPdxi(1,ijn)+fxp*dPdxi(1,ijp))*xpn*nxx
  ! dpyii = (fxn*dPdxi(2,ijn)+fxp*dPdxi(2,ijp))*ypn*nyy
  ! dpzii = (fxn*dPdxi(3,ijn)+fxp*dPdxi(3,ijp))*zpn*nzz
  ! flmass = dene*(ui*arx+vi*ary+wi*arz) + cap*(p(ijn)-p(ijp)-dpxii-dpyii-dpzii)

end subroutine
