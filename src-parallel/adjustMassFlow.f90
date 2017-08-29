subroutine adjustMassFlow
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables

  implicit none

  integer :: i,ijb,ijp,iface
  real(dp) :: flowo,fac


  ! Inlet boundaries (mass fluxes prescribed in routine 'bcin')
  do i=1,ninl

    iface = iInletFacesStart+i
    ijp = owner(iface)

    ! Minus sign is there to make fmi(i) positive since it enters the cell.
    ! Check out comments in bcin.f90
    su(ijp) = su(ijp) - fmi(i) 

  end do


  ! Loop Outlet cells, first to get flowo, then again after calculating fac.

  ! Extrapolated velocity at outlet boundary, outlet mass fluxes
  flowo=0.0_dp

  do i=1,nout

    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i

    u(ijb) = u(ijp)
    v(ijb) = v(ijp)
    w(ijb) = w(ijp)

    fmo(i) = den(ijp)*( u(ijb)*arx(iface)+v(ijb)*ary(iface)+w(ijb)*arz(iface) )
    

    flowo = flowo + fmo(i) 

  end do

  ! Correct mass flux to satisfy global mass conservation & add to source
  fac = flomas/(flowo+small)

  do i=1,nout
    
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i


    fmo(i) = fmo(i)*fac

    u(ijb) = u(ijb)*fac
    v(ijb) = v(ijb)*fac
    w(ijb) = w(ijb)*fac

    ! fmo is positive because of same direction of velocity and surface normal vectors
    ! but the mass flow is going out of the cell, therefore minus sign.
    su(ijp) = su(ijp) - fmo(i)

  end do

end subroutine