subroutine correctBoundaryConditionsVelocity
!  
!******************************************************************************
!
!     Updates values at symmetry boundaries 
! 
!******************************************************************************
!
  use types
  use parameters
  use geometry
  use variables

  implicit none

!
!     Local variables
!
  integer :: i,ijp,ijb,iface
  real(dp) :: Unmag


  ! Update velocity components along symmetry boundaries
  do i=1,nsym

    iface = iSymmetryFacesStart+i
    ijp = owner(iface)
    ijb = iSymmetryStart+i

    ! Project velocity vector to face normal direction:
    Unmag = u(ijp)*arx(iface)+v(ijp)*ary(iface)+w(ijp)*arz(iface)

    U(ijb) = U(ijp)-Unmag*arx(iface)
    V(ijb) = V(ijp)-Unmag*ary(iface)
    W(ijb) = W(ijp)-Unmag*arz(iface)

  end do

end subroutine
