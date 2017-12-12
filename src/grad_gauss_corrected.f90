subroutine grad_gauss_corrected(u,dudx,dudy,dudz)
!
!***********************************************************************
!
!     Calculates cell centered gradient using gauss theorem
!     parameters
!     u - field, the gradient of which we are looking for
!     dudx,dudy,dudz - arrays where the gradient components are stored
!
!     gauss gradient rule:
!     ------->                                 ->
!     grad(u) = 1/vol * sum_{i=1}^{i=nf} (u)_f*sf
!     where:
!     grad(u) - cell centered gradient vector
!     (u)_f   - face interpolated value of scalar u
!     vol     - cell volume
!     sf      - cell face area vector
!     nf      - number of faces in a cell
!
!***********************************************************************
!
  use types
  use parameters
  use geometry

  implicit none

  ! Arguments
  real(dp), dimension(numTotal), intent(in) :: u
  real(dp), dimension(numCells), intent(inout) :: dudx,dudy,dudz

  ! Local
  integer :: i,ijp,ijn,ijb,iface
  real(dp) :: volr
  real(dp), dimension(numCells) :: dfxo,dfyo,dfzo

  ! Initialize gradient with lsq gradient
  dfxo = dudx
  dfyo = dudy
  dfzo = dudz
  
  ! Initialize new gradient
  dudx = 0.0_dp
  dudy = 0.0_dp
  dudz = 0.0_dp

  ! Calculate terms integrated over surfaces

  ! Inner face
  do i=1,numInnerFaces
    ijp = owner(i)
    ijn = neighbour(i)
    call gradco(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                u, dfxo, dfyo, dfzo, dudx, dudy, dudz)
  enddo

  ! Contribution from O- and C-grid cuts
  do i=1,noc
    iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
    ijp = ijl(i)
    ijn = ijr(i)
    call gradco(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), &
                u, dfxo, dfyo, dfzo, dudx, dudy, dudz)
  end do

  ! Contribution from boundaries
  do i = 1,ninl
    iface = iInletFacesStart+i
    ijp = owner(iface)
    ijb = iInletStart + i
    call gradbc(arx(iface), ary(iface), arz(iface), u(ijb), dudx(ijp), dudy(ijp), dudz(ijp))
  enddo

  do i = 1,nout
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart + i
    call gradbc(arx(iface), ary(iface), arz(iface), u(ijb), dudx(ijp), dudy(ijp), dudz(ijp))
  enddo

  do i = 1,nsym
    iface = iSymmetryFacesStart+i
    ijp = owner(iface)
    ijb = iSymmetryStart+i
    call gradbc(arx(iface), ary(iface), arz(iface), u(ijb), dudx(ijp), dudy(ijp), dudz(ijp))
  enddo   

  do i = 1,nwal
    iface = iWallFacesStart+i
    ijp = owner(iface)
    ijb = iWallStart+i
    call gradbc(arx(iface), ary(iface), arz(iface), u(ijb), dudx(ijp), dudy(ijp), dudz(ijp))
  enddo

  do i=1,npru
    iface = iPressOutletFacesStart + i
    ijp = owner(iface)
    ijb = iPressOutletStart + i
    call gradbc(arx(iface), ary(iface), arz(iface), u(ijb), dudx(ijp), dudy(ijp), dudz(ijp))
  enddo


  ! Calculate gradient components at cv-centers
  do ijp=1,numCells
    volr=1.0_dp/vol(ijp)
    dudx(ijp)=dudx(ijp)*volr
    dudy(ijp)=dudy(ijp)*volr
    dudz(ijp)=dudz(ijp)*volr
  enddo

return
end