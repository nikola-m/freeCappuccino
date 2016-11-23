!***********************************************************************
!
subroutine Additional_algebraic_heatflux_terms
!
!***********************************************************************
!     Calculates the additional agebraic heatflux terms for temperature eq.
!     and ads them to SU rhs. vector.
!
!***********************************************************************
  use types, only: dp
  use parameters, only: numInnerFaces,viscos
  use indexes, only: owner,neighbour
  use geometry, only: facint,arx,ary,arz
  use sparse_matrix, only: su
  use variables, only: vis,den
  use temperature, only: dTdxi, utt,vtt,wtt

  implicit none
!
!***********************************************************************
!

!
!     Local variables
!
  integer :: i, ijp, ijn
  real(dp) :: term1i, term2i,term3i
  real(dp) :: fxp, fxn
  real(dp) :: prant,prtr

  prant=0.86_dp
  prtr=1./prant

!-------------------------------------------------------
!
!     TERM1 is interpolated to faces, then a x-direction
!     gradient is sought for.
!     TERM1 is \rho*utt + dTdx*prtr*\mu_t
!     where: 
!       \mu_t = \mu_eff - \mu_molecular,
!       utt - heat flux
!
!     TERM2 is interpolated to faces, then a y-direction
!     gradient is sought for.
!     TERM2 is \rho*vtt + dTdy*prtr*\mu_t
!
!     TERM3 is interpolated to faces, then a z-direction
!     gradient is sought for.
!     TERM3 is \rho*wtt + dTdz*prtr*\mu_t
!
!     su(inp)=su(inp)-dterm1dx-dterm2dy-dterm3dz
!-------------------------------------------------------

do i=1,numInnerFaces
  ijp = owner(i)
  ijn = neighbour(i)

  fxn = facint(ijp)
  fxp = 1.0_dp-facint(ijp)

        ! Interpolate term1:
        term1i = (den(ijn)*utt(ijn)+dTdxi(1,ijn)*prtr*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*utt(ijp)+dTdxi(1,ijp)*prtr*(vis(ijp)-viscos))*fxp

        ! Interpolate term2:
        term2i = (den(ijn)*vtt(ijn)+dTdxi(2,ijn)*prtr*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*vtt(ijp)+dTdxi(2,ijp)*prtr*(vis(ijp)-viscos))*fxp

        ! Interpolate term3:
        term3i = (den(ijn)*wtt(ijn)+dTdxi(3,ijn)*prtr*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*wtt(ijp)+dTdxi(3,ijp)*prtr*(vis(ijp)-viscos))*fxp

        su(ijp) = su(ijp) - ( term1i * arx(i) + term2i * ary(i) + term3i * arz(i) )

        su(ijn) = su(ijn) + ( term1i * arx(i) + term2i * ary(i) + term3i * arz(i) )
enddo


end subroutine