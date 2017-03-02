!***********************************************************************
!
subroutine fluxmc(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, fmcor)
!
!***********************************************************************
!
!   This routine calculates mass flux correction in the
!   second pressure-correction step which accounts for the
!   effects of non-orthogonality
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: xc,yc,zc,vol
  use variables
  use sparse_matrix, only: apu

  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), intent(inout) :: fmcor
  !
  ! Local variables
  !
  real(dp) :: fxn,fxp
  real(dp) :: xpn,ypn,zpn
  real(dp) :: rapr
  real(dp) :: are
  real(dp) :: nxx,nyy,nzz
  real(dp) :: xep,yep,zep,xpp,ypp,zpp
  real(dp) :: dppnnr

  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

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

  ! Distance from P' to N'-reciprocal value
  dppnnr = 1.0_dp/((xpn*nxx)+(ypn*nyy)+(zpn*nzz))

  ! Values at points p' and e' due to non-orthogonality. 
  xpp=xf-(xf-xc(ijp))*nxx
  ypp=yf-(yf-yc(ijp))*nyy
  zpp=zf-(zf-zc(ijp))*nzz

  xep=xf-(xf-xc(ijn))*nxx
  yep=yf-(yf-yc(ijn))*nyy
  zep=zf-(zf-zc(ijn))*nzz

  ! Distances |P'P| and |E'E| projected ionto x,y,z-axis
  xpp=xpp-xc(ijp)
  ypp=ypp-yc(ijp)
  zpp=zpp-zc(ijp)

  xep=xep-xc(ijn)
  yep=yep-yc(ijn)
  zep=zep-zc(ijn)

  ! APU==1./AP x density - interpolated            
  rapr = (apu(ijp)*den(ijp)*vol(ijp)*fxp+apu(ijn)*den(ijn)*vol(ijn)*fxn)


  ! Mass flux correction for the second p'-equation (source term)
  fmcor = rapr*are*((dPdxi(1,ijn)*xep-dPdxi(1,ijp)*xpp)   & 
                   +(dPdxi(2,ijn)*yep-dPdxi(2,ijp)*ypp)   &
                   +(dPdxi(3,ijn)*zep-dPdxi(3,ijp)*zpp))  &
                   *dppnnr 

end subroutine
