subroutine adjustMassFlow
  use types
  use parameters
  use indexes
  use geometry
  use sparse_matrix
  use variables

  implicit none

  integer :: i,ijb,ijp
  real(dp) :: flowo,fac


  ! Inlet boundaries (mass fluxes prescribed in routine 'bcin')
  do i=1,ninl

    ijp = owner(iInletFacesStart+i)

    ! Minus sign is there to make fmi(i) positive since it enters the cell.
    ! Check out comments in bcin.f90
    su(ijp) = su(ijp) - fmi(i) 

  end do


  ! Loop Outlet cells, first to get flowo, then again after calculating fac.

  ! Extrapolated velocity at outlet boundary, outlet mass fluxes
  flowo=0.0d0

  do i=1,nout

    ijp = owner(iOutletFacesStart+i)
    ijb = iOutletStart+i

    u(ijb) = u(ijp)
    v(ijb) = v(ijp)
    w(ijb) = w(ijp)

    fmo(i) = den(ijp)*( u(ijb)*xno(i)+v(ijb)*yno(i)+w(ijb)*zno(i) )
    

    flowo = flowo + fmo(i) 

  end do

  ! Correct mass flux to satisfy global mass conservation & add to source
  fac = flomas/(flowo+small)

  do i=1,nout

    ijp = owner(iOutletFacesStart+i)
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