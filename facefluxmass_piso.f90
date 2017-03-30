!***********************************************************************
!
subroutine facefluxmass_piso(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, cap, can, flmass)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: xc,yc,zc,vol
  use sparse_matrix, only: apu
  use variables, only: den,U,V,W,dUdxi,dVdxi,dWdxi
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
  real(dp) :: are,dpn
  real(dp) :: xpn,ypn,zpn,dene,smdpn
  real(dp) :: xi,yi,zi
  real(dp) :: nxx,nyy,nzz
  real(dp) :: ui,vi,wi
  real(dp) :: duxi,duyi,duzi


  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Coordinates of point e'
  xi=xc(ijp)*fxp+xc(ijn)*fxn
  yi=yc(ijp)*fxp+yc(ijn)*fxn
  zi=zc(ijp)*fxp+zc(ijn)*fxn

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2) 

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! Unit vectors of the normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are


  ! density at the cell face
  dene=den(ijp)*fxp+den(ijn)*fxn

  ! COEFFICIENTS OF PRESSURE EQUATION
  ! sfdpnr=1./(ARX*XPN*nxx+ARY*YPN*nyy+ARZ*ZPN*nzz)
  ! smdpn = (arx*arx+ary*ary+arz*arz)/(arx*xpn+ary*ypn+arz*zpn)
  smdpn = are/dpn
  cap = -dene*(fxp*vol(ijp)*apu(ijp)+fxn*vol(ijn)*apu(ijn))*smdpn
  can = cap


  !+Interpolate velocities to face center:+++++++++++++++++++++++++
  !    |________Ue'_________|
  ui = u(ijp)*fxp+u(ijn)*fxn

  ! Interpolate gradients defined at CV centers to faces
  duxi = dudxi(1,ijp)*fxp+dudxi(1,ijn)*fxn
  duyi = dudxi(2,ijp)*fxp+dudxi(2,ijn)*fxn
  duzi = dudxi(3,ijp)*fxp+dudxi(3,ijn)*fxn
  !      |_______________Ucorr________________|
  ui = ui+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))



  !  |________Ve'_________|
  vi=v(ijp)*fxp+v(ijn)*fxn

  duxi = dvdxi(1,ijp)*fxp+dvdxi(1,ijn)*fxn
  duyi = dvdxi(2,ijp)*fxp+dvdxi(2,ijn)*fxn
  duzi = dvdxi(3,ijp)*fxp+dvdxi(3,ijn)*fxn
  !    |_______________Vcorr________________|
  vi=vi+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi))

  

 !  |________We'_________|
  wi=w(ijp)*fxp+w(ijn)*fxn

  duxi = dwdxi(1,ijp)*fxp+dwdxi(1,ijn)*fxn
  duyi = dwdxi(2,ijp)*fxp+dwdxi(2,ijn)*fxn
  duzi = dwdxi(3,ijp)*fxp+dwdxi(3,ijn)*fxn
  !    |_______________Wcorr________________|
  wi=wi+(duxi*(xf-xi)+duyi*(yf-yi)+duzi*(zf-zi)) 
  
  !+END: Interpolate velocities to face center:+++++++++++++++++++++++++

 
  ! MASS FLUX
  !// calculate the fluxes by dotting the interpolated velocity (to cell faces) with face normals
  !     phi = (fvc::interpolate(U) & mesh.Sf()) 
  flmass=dene*(ui*arx+vi*ary+wi*arz)

end subroutine
