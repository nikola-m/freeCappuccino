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
  integer :: i,ijp,ijb
  real(dp) :: Unmag


  ! Update velocity components along symmetry boundaries
  do i=1,nsym

    ijp = owner(iSymmetryFacesStart+i)
    ijb = iSymmetryStart+i

    ! Project velocity vector to face normal direction:
    Unmag = u(ijp)*Xns(i)+v(ijp)*Yns(i)+w(ijp)*Zns(i)

    U(ijb) = U(ijp)-Unmag*xns(i)
    V(ijb) = V(ijp)-Unmag*yns(i)
    W(ijb) = W(ijp)-Unmag*zns(i)

  end do

end subroutine
