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
  use interpolation

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
  real(dp) :: nxx,nyy,nzz
  real(dp) :: ui,vi,wi


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


  ! Interpolate velocities to face center:

  ! ui = face_value_cds_corrected( ijp, ijn, xf, yf, zf, lambda, u, dUdxi )
  ! vi = face_value_cds_corrected( ijp, ijn, xf, yf, zf, lambda, v, dVdxi )
  ! wi = face_value_cds_corrected( ijp, ijn, xf, yf, zf, lambda, w, dWdxi )


  ui = face_value_central( ijp,ijn, xf, yf, zf, u, dUdxi )
  vi = face_value_central( ijp,ijn, xf, yf, zf, v, dVdxi )
  wi = face_value_central( ijp,ijn, xf, yf, zf, w, dWdxi )


  ! MASS FLUX
  !// calculate the fluxes by dotting the interpolated velocity (to cell faces) with face normals
  !     phi = (fvc::interpolate(U) & mesh.Sf()) 
  flmass = dene*(ui*arx+vi*ary+wi*arz)

end subroutine
