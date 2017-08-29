!***********************************************************************
!
subroutine calcstress
!***********************************************************************
!
! Calculates turbulent stresses
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: numCells, vol
  use variables

  implicit none 
!
!***********************************************************************
!
  integer :: inp
  real(dp) :: volr                    ! 1/cellvolume
  real(dp) :: vist                    ! turbulent viscosity
  real(dp) :: dudx, dudy, dudz, &     ! dudxi - the velocity gradient
              dvdx, dvdy, dvdz, &     ! dvdxi - the velocity gradient
              dwdx, dwdy, dwdz        ! dwdxi - the velocity gradient
  real(dp) ::  uuold,vvold, wwold, &  ! Reynolds stress tensor  components 
              uvold, uwold,vwold 
  real(dp) :: facnapm

  facnapm = 1.0_dp-facnap

  do inp=1,numCells

    volr=1./vol(inp)

    vist=(vis(inp)-viscos)/densit

    uuold=uu(inp)
    vvold=vv(inp)
    wwold=ww(inp)
    uvold=uv(inp)
    uwold=uw(inp)
    vwold=vw(inp)

    dudx = dUdxi(1,inp)
    dudy = dUdxi(2,inp)
    dudz = dUdxi(3,inp)
    
    dvdx = dVdxi(1,inp)
    dvdy = dVdxi(2,inp)
    dvdz = dVdxi(3,inp)

    dwdx = dWdxi(1,inp)
    dwdy = dWdxi(2,inp)
    dwdz = dWdxi(3,inp)
!
!==============================================
!---[EVM approach: ]
!==============================================
    if(levm) then

      uu(inp)=twothirds*te(inp)-vist*(dudx+dudx)
      vv(inp)=twothirds*te(inp)-vist*(dvdy+dvdy)
      ww(inp)=twothirds*te(inp)-vist*(dwdz+dwdz)

      uv(inp)=-vist*(dudy+dvdx)
      uw(inp)=-vist*(dudz+dwdx)
      vw(inp)=-vist*(dvdz+dwdy)
!
!==============================================

!==============================================
!-----------------------------------------------------
!     don't forget to include exact production GEN !!
!-----------------------------------------------------
    else if(lasm) then

      ! Wallin-Johansson EARSM
      uu(inp) = twothirds*te(inp)-vist*(dudx+dudx) + 2.*bij(1,inp)*te(inp)
      vv(inp) = twothirds*te(inp)-vist*(dvdy+dvdy) + 2.*bij(4,inp)*te(inp)
      ! It seems that b(3,3)=-b(1,1)-b(2,2):
      ww(inp) = twothirds*te(inp)-vist*(dwdz+dwdz)  &
              - 2.*(bij(1,inp)+bij(4,inp))*te(inp) 

      uv(inp) = -vist*(dudy+dvdx) + 2.*bij(2,inp)*te(inp)
      uw(inp) = -vist*(dudz+dwdx) + 2.*bij(3,inp)*te(inp)
      vw(inp) = -vist*(dvdz+dwdy) + 2.*bij(5,inp)*te(inp)
 
! !==========================================================
! !     Reynolds stresses using the anisotropy tensor:
! !     ______
! !     u_iu_j = k(a_ij + 2/3 \delta_ij)
! !     a_ij = 2* b_ij   
! !==========================================================
!       uu(inp)=2*bij(1,inp)*te(inp) + 2*te(inp)/3.
!       vv(inp)=2*bij(4,inp)*te(inp) + 2*te(inp)/3.
!       ww(inp)=0.                   + 2*te(inp)/3.

!       uv(inp)=2*bij(2,inp)*te(inp)
!       uw(inp)=2*bij(3,inp)*te(inp)
!       vw(inp)=2*bij(5,inp)*te(inp)

    end if !![for evm or asm approach]
!----------------------------------------------

    ! Clip negative values
    uu(inp)=max(uu(inp),small)
    vv(inp)=max(vv(inp),small)
    ww(inp)=max(ww(inp),small)

    ! Underrelax using facnap factor
    uu(inp)=facnap*uu(inp)+facnapm*uuold
    vv(inp)=facnap*vv(inp)+facnapm*vvold
    ww(inp)=facnap*ww(inp)+facnapm*wwold

    uv(inp)=facnap*uv(inp)+facnapm*uvold
    uw(inp)=facnap*uw(inp)+facnapm*uwold
    vw(inp)=facnap*vw(inp)+facnapm*vwold


  end do ! cell-loop
!-------------------------------

end subroutine
