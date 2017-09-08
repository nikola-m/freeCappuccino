!***********************************************************************
!
subroutine Additional_algebraic_stress_terms
!
!***********************************************************************
!     Calculates the additional agebraic stress terms for momentum eq.
!     and ads them to SU, SV, SW rhs. vectors.
!
!***********************************************************************
  use types, only: dp
  use parameters, only: viscos
  use geometry, only: numCells,numInnerFaces,npro,facint,fpro,arx,ary,arz,owner,neighbour,iProcFacesStart,iProcStart
  use sparse_matrix, only: su,sv,sw
  use variables, only: den,vis,uu,uv,uw,vv,vw,ww,dudxi,dvdxi,dwdxi,dtedxi

  implicit none
!
!***********************************************************************
!

!
!     Local variables
!
  integer :: i, ijp, ijn, inp, iface
  real(dp) :: term1i, term2i,term3i
  real(dp) :: two_thirds, fxp, fxn


  two_thirds=2./3.0_dp

!-------------------------------------------------------
!     [U-VELOCITY COMPONENT: ]
!
!     TERM1 is interpolated to faces, then a x-direction
!     gradient is sought for.
!     TERM1 is \rho*tau_{1,1} + 2*dUdx*\mu_t
!     where: 
!       \mu_t = \mu_eff - \mu_molecular,
!       tau_{1,i} - first row of Reynolds stress tensor
!
!     TERM2 is interpolated to faces, then a y-direction
!     gradient is sought for.
!     TERM2 is \rho*tau_{1,2} + \mu_t*(dUdy+DVdx)
!
!     TERM3 is interpolated to faces, then a z-direction
!     gradient is sought for.
!     TERM3 is \rho*tau_{1,3} + \mu_t*(dUdz+DWdx)
!
!     su(inp)=su(inp)-dterm1dx-dterm2dy-dterm3dz+ two_thirds*dtedxi(1,inp)
!-------------------------------------------------------

do i=1,numInnerFaces
  ijp = owner(i)
  ijn = neighbour(i)

  fxn = facint(ijp)
  fxp = 1.0_dp-fxn

        ! Interpolate term1:
        term1i = (den(ijn)*uu(ijn)+(dudxi(1,ijn)+dudxi(1,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*uu(ijp)+(dudxi(1,ijp)+dudxi(1,ijp))*(vis(ijp)-viscos))*fxp

        ! Interpolate term2:
        term2i = (den(ijn)*uv(ijn)+(dudxi(2,ijn)+dvdxi(1,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*uv(ijp)+(dudxi(2,ijp)+dvdxi(1,ijp))*(vis(ijp)-viscos))*fxp

        ! Interpolate term3:
        term3i = (den(ijn)*uw(ijn)+(dudxi(3,ijn)+dwdxi(1,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*uw(ijp)+(dudxi(3,ijp)+dwdxi(1,ijp))*(vis(ijp)-viscos))*fxp

        su(ijp) = su(ijp) - ( term1i * arx(i) + term2i * ary(i) + term3i * arz(i) )

        su(ijn) = su(ijn) + ( term1i * arx(i) + term2i * ary(i) + term3i * arz(i) )
enddo

! Faces at process boundaries
do i=1,npro
  iface = iProcFacesStart + i
  ijp = owner(iface)
  ijn = iProcStart + i

  fxn = fpro(ijp)
  fxp = 1.0_dp-fxn

        ! Interpolate term1:
        term1i = (den(ijn)*uu(ijn)+(dudxi(1,ijn)+dudxi(1,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*uu(ijp)+(dudxi(1,ijp)+dudxi(1,ijp))*(vis(ijp)-viscos))*fxp

        ! Interpolate term2:
        term2i = (den(ijn)*uv(ijn)+(dudxi(2,ijn)+dvdxi(1,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*uv(ijp)+(dudxi(2,ijp)+dvdxi(1,ijp))*(vis(ijp)-viscos))*fxp

        ! Interpolate term3:
        term3i = (den(ijn)*uw(ijn)+(dudxi(3,ijn)+dwdxi(1,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*uw(ijp)+(dudxi(3,ijp)+dwdxi(1,ijp))*(vis(ijp)-viscos))*fxp

        su(ijp) = su(ijp) - ( term1i * arx(iface) + term2i * ary(iface) + term3i * arz(iface) )
enddo

do inp=1,numCells
   su(inp) = su(inp) + two_thirds*dtedxi(1,inp)
end do


!...sada isto za V komponentu
do i=1,numInnerFaces
  ijp = owner(i)
  ijn = neighbour(i)

  fxn = facint(ijp)
  fxp = 1.0_dp-fxn

        ! Interpolate term1:
        term1i = (den(ijn)*uv(ijn)+(dvdxi(1,ijn)+dudxi(2,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*uv(ijp)+(dvdxi(1,ijp)+dudxi(2,ijp))*(vis(ijp)-viscos))*fxp

        ! Interpolate term2:
        term2i = (den(ijn)*vv(ijn)+(dvdxi(2,ijn)+dvdxi(2,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*vv(ijp)+(dvdxi(2,ijp)+dvdxi(2,ijp))*(vis(ijp)-viscos))*fxp

        ! Interpolate term3:
        term3i = (den(ijn)*vw(ijn)+(dvdxi(3,ijn)+dwdxi(2,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*vw(ijp)+(dudxi(3,ijp)+dwdxi(1,ijp))*(vis(ijp)-viscos))*fxp

        sv(ijp) = sv(ijp) - ( term1i * arx(i) + term2i * ary(i) + term3i * arz(i) )

        sv(ijn) = sv(ijn) + ( term1i * arx(i) + term2i * ary(i) + term3i * arz(i) )
enddo


! Faces at process boundaries
do i=1,npro
  iface = iProcFacesStart + i
  ijp = owner(iface)
  ijn = iProcStart + i

  fxn = fpro(ijp)
  fxp = 1.0_dp-fxn

        ! Interpolate term1:
        term1i = (den(ijn)*uv(ijn)+(dvdxi(1,ijn)+dudxi(2,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*uv(ijp)+(dvdxi(1,ijp)+dudxi(2,ijp))*(vis(ijp)-viscos))*fxp

        ! Interpolate term2:
        term2i = (den(ijn)*vv(ijn)+(dvdxi(2,ijn)+dvdxi(2,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*vv(ijp)+(dvdxi(2,ijp)+dvdxi(2,ijp))*(vis(ijp)-viscos))*fxp

        ! Interpolate term3:
        term3i = (den(ijn)*vw(ijn)+(dvdxi(3,ijn)+dwdxi(2,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*vw(ijp)+(dudxi(3,ijp)+dwdxi(1,ijp))*(vis(ijp)-viscos))*fxp

        sv(ijp) = sv(ijp) - ( term1i * arx(iface) + term2i * ary(iface) + term3i * arz(iface) )

enddo


do inp=1,numCells
   sv(inp) = sv(inp) + two_thirds*dtedxi(2,inp)
end do


!...sada isto za W komponentu
do i=1,numInnerFaces
  ijp = owner(i)
  ijn = neighbour(i)

  fxn = facint(ijp)
  fxp = 1.0_dp-fxn

        ! Interpolate term1:
        term1i = (den(ijn)*uw(ijn)+(dwdxi(1,ijn)+dudxi(3,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(inp)*uw(inp)+(dwdxi(1,inp)+dudxi(3,inp))*(vis(inp)-viscos))*fxp

        ! Interpolate term2:
        term2i = (den(ijn)*vw(ijn)+(dwdxi(2,ijn)+dvdxi(3,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*vw(ijp)+(dwdxi(2,ijp)+dvdxi(3,ijp))*(vis(ijp)-viscos))*fxp

        ! Interpolate term3:
        term3i = (den(ijn)*ww(ijn)+(dwdxi(3,ijn)+dwdxi(3,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*ww(ijp)+(dwdxi(3,ijp)+dwdxi(3,ijp))*(vis(ijp)-viscos))*fxp

        sw(ijp) = sw(ijp) - ( term1i * arx(i) + term2i * ary(i) + term3i * arz(i) )

        sw(ijn) = sw(ijn) + ( term1i * arx(i) + term2i * ary(i) + term3i * arz(i) )
enddo

! Faces at process boundaries
do i=1,npro
  iface = iProcFacesStart + i
  ijp = owner(iface)
  ijn = iProcStart + i

  fxn = fpro(ijp)
  fxp = 1.0_dp-fxn

        ! Interpolate term1:
        term1i = (den(ijn)*uw(ijn)+(dwdxi(1,ijn)+dudxi(3,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(inp)*uw(inp)+(dwdxi(1,inp)+dudxi(3,inp))*(vis(inp)-viscos))*fxp

        ! Interpolate term2:
        term2i = (den(ijn)*vw(ijn)+(dwdxi(2,ijn)+dvdxi(3,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*vw(ijp)+(dwdxi(2,ijp)+dvdxi(3,ijp))*(vis(ijp)-viscos))*fxp

        ! Interpolate term3:
        term3i = (den(ijn)*ww(ijn)+(dwdxi(3,ijn)+dwdxi(3,ijn))*(vis(ijn)-viscos))*fxn+ &
                 (den(ijp)*ww(ijp)+(dwdxi(3,ijp)+dwdxi(3,ijp))*(vis(ijp)-viscos))*fxp

        sw(ijp) = sw(ijp) - ( term1i * arx(iface) + term2i * ary(iface) + term3i * arz(iface) )

enddo

do inp=1,numCells
   sw(inp) = sw(inp) + two_thirds*dtedxi(3,inp)
end do


end subroutine
