subroutine grad_gauss(u,dudx,dudy,dudz)
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

  !**Arguments************************************************************
  real(dp), dimension(numTotal), intent(in) :: u
  real(dp), dimension(numCells), intent(inout) :: dudx,dudy,dudz

  !**Locals***************************************************************
  integer :: i,ijp,ijn,ijb,lc,iface
  real(dp) :: volr
  real(dp), dimension(numCells) :: dfxo,dfyo,dfzo

  ! Initialize gradient
  dfxo(:) = 0.0d0
  dfyo(:) = 0.0d0
  dfzo(:) = 0.0d0

  ! Start iterative calculation of gradients
  do lc = 1,nigrad
    
    ! Initialize new gradient
    dudx(:) = 0.0d0
    dudy(:) = 0.0d0
    dudz(:) = 0.0d0

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
      iface = iOCFacesStart+i
      ijp = ijl(i)
      ijn = ijr(i)
      call gradco(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), &
                  u, dfxo, dfyo, dfzo, dudx, dudy, dudz)
    end do

    ! Contribution from boundaries
    do i=1,numBoundaryFaces
      iface = numInnerFaces + i
      ijp = owner(iface)
      ijb = numCells+i
      call gradbc(ijp, ijb, arx(iface), ary(iface), arz(iface), u, dudx, dudy, dudz)
    end do 

    ! Calculate gradient components at cv-centers
    do ijp=1,numCells
          volr=1./vol(ijp)
          dudx(ijp)=dudx(ijp)*volr
          dudy(ijp)=dudy(ijp)*volr
          dudz(ijp)=dudz(ijp)*volr
          if(lc.eq.nigrad) print*,ijp,dudx(ijp)
    enddo

    ! Set old gradient = new gradient for the next iteration
    if(lc.ne.nigrad) then
      dfxo=dudx
      dfyo=dudy
      dfzo=dudz
    endif

  enddo ! lc-loop

end subroutine



subroutine gradco(ijp,ijn, &
                  xfc,yfc,zfc,sx,sy,sz,fif, &
                  fi,dfxo,dfyo,dfzo,dfx,dfy,dfz)
!=======================================================================
!     This routine calculates contribution to the gradient
!     vector of a scalar FI at the CV center, arising from
!     an inner cell face (cell-face value of FI times the 
!     corresponding component of the surface vector).
!=======================================================================
  use types
  use parameters
  use geometry

  implicit none

  integer,    intent(in) :: ijp,ijn
  real(dp), intent(in) :: xfc,yfc,zfc
  real(dp), intent(in) :: sx,sy,sz
  real(dp), intent(in) :: fif
  real(dp), dimension(numTotal), intent(in) :: fi
  real(dp), dimension(numCells), intent(in) :: dfxo,dfyo,dfzo
  real(dp), dimension(numCells), intent(inout)  :: dfx,dfy,dfz


  real(dp) :: xi,yi,zi,dfxi,dfyi,dfzi
  real(dp) :: fie,dfxe,dfye,dfze
  real(dp) :: fxn,fxp

    !
    ! Coordinates of point on the line connecting center and neighbor,
    ! old gradient vector components interpolated for this location.

    fxn=fif 
    fxp=1.0d0-fxn

    xi=xc(ijp)*fxp+xc(ijn)*fxn
    yi=yc(ijp)*fxp+yc(ijn)*fxn
    zi=zc(ijp)*fxp+zc(ijn)*fxn

    dfxi=dfxo(ijp)*fxp+dfxo(ijn)*fxn
    dfyi=dfyo(ijp)*fxp+dfyo(ijn)*fxn
    dfzi=dfzo(ijp)*fxp+dfzo(ijn)*fxn

    ! Value of the variable at cell-face center
    fie=fi(ijp)*fxp+fi(ijn)*fxn + dfxi*(xfc-xi)+dfyi*(yfc-yi)+dfzi*(zfc-zi)

    ! (interpolated mid-face value)x(area)
    dfxe=fie*sx
    dfye=fie*sy
    dfze=fie*sz
    
    ! Accumulate contribution at cell center and neighbor
    dfx(ijp)=dfx(ijp)+dfxe
    dfy(ijp)=dfy(ijp)+dfye
    dfz(ijp)=dfz(ijp)+dfze
    !--
    dfx(ijn)=dfx(ijn)-dfxe
    dfy(ijn)=dfy(ijn)-dfye
    dfz(ijn)=dfz(ijn)-dfze

end subroutine



subroutine gradbc(ijp,ijb,sx,sy,sz,fi,dfx,dfy,dfz)
!=======================================================================
!     This routine calculates the contribution of a 
!     boundary cell face to the gradient at CV-center.
!=======================================================================
  use types
  use geometry, only: numTotal,numCells

  implicit none

  integer,    intent(in) :: ijp,ijb
  real(dp), intent(in) :: sx,sy,sz
  real(dp), dimension(numTotal), intent(in) :: fi
  real(dp), dimension(numCells), intent(inout)  :: dfx,dfy,dfz

  dfx(ijp)=dfx(ijp)+fi(ijb)*sx
  dfy(ijp)=dfy(ijp)+fi(ijb)*sy
  dfz(ijp)=dfz(ijp)+fi(ijb)*sz

end subroutine