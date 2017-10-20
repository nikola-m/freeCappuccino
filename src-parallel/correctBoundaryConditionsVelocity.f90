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
  real(dp) :: Unmag,flowo,fac

  ! Update velocity components along outlet boundaries
  ! and correct mass flux to satisfy global mass conservation

  flowo=0.0_dp

  do i=1,nout

    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i

    U(ijb) = U(ijp)
    V(ijb) = V(ijp)
    W(ijb) = W(ijp)

    fmo(i) = den(ijp)*( u(ijb)*arx(iface)+v(ijb)*ary(iface)+w(ijb)*arz(iface) )
    
    flowo = flowo + fmo(i)

  end do

  ! Ratio of inflow and outflow mass flux
  fac = flomas/(flowo+small)

  ! Correct mass flux to satisfy global mass conservation
  do i=1,nout
    
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i


    fmo(i) = fmo(i)*fac

    u(ijb) = u(ijb)*fac
    v(ijb) = v(ijb)*fac
    w(ijb) = w(ijb)*fac

  end do

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
