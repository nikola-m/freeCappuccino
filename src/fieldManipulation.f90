    module fieldManipulation
!
!***********************************************************************
!
  use types
  use parameters
  use gradients

  implicit none

  public

  contains

!***********************************************************************
!
  pure function volumeWeightedAverage(U) result(wAvgU)
!
!***********************************************************************
    use geometry, only:numTotal,numCells,vol

    implicit none

!...Output
    real(dp) :: wAvgU

!...Input
    real(dp), dimension(numTotal), intent(in) :: U

!...Locals
    integer :: inp
    real(dp) :: sumvol
!
!***********************************************************************
!  
    sumvol = 0.0_dp
    wAvgU = 0.0_dp 

      do inp=1,numCells
          wAvgU = wAvgU + (Vol(inp)*U(inp))
          sumvol = sumvol + vol(inp)
      enddo
    
    wAvgU = wAvgU / sumvol

  end function


!***********************************************************************
!
  subroutine calcPressDiv 
!
!***********************************************************************
!  
!  -fvc:Div(p)                                       
!  ExplDiv(u) = sum_{i=1}^{i=nf} (u)_f*sf or
!  Interpolation to cell face centers done by cds corrected scheme.
!
!***********************************************************************
!
  use geometry
  use variables, only: p,dPdxi
  use sparse_matrix, only: su,sv,sw

  implicit none
!
!***********************************************************************
!


!...Local
    integer :: i,ijp,ijn,ijb,iface,istage
    real(dp) :: dfxe,dfye,dfze

    ! Pressure gradient
    do istage=1,nipgrad
      ! Pressure at boundaries (for correct calculation of press. gradient)
      call bpres(p,istage)
      ! Calculate pressure gradient.
      call grad(p,dPdxi)
    end do

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call presFaceDivInner(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                        p, dPdxi,dfxe,dfye,dfze)

      ! Accumulate contribution at cell center and neighbour.
      ! Note, we calculate negative Dievrgence, therefore opposite sign...
      su(ijp) = su(ijp)-dfxe
      sv(ijp) = sv(ijp)-dfye
      sw(ijp) = sw(ijp)-dfze
       
      su(ijn) = su(ijn)+dfxe
      sv(ijn) = sv(ijn)+dfye
      sw(ijn) = sw(ijn)+dfze

    enddo

    ! Contribution from O- and C-grid cuts
    do i=1,noc
      iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
      ijp = ijl(i)
      ijn = ijr(i)
      call presFaceDivInner(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), &
                        p, dPdxi,dfxe,dfye,dfze)

      ! Accumulate contribution at cell center and neighbour.
      ! Note, we calculate negative Dievrgence, therefore opposite sign...
      su(ijp) = su(ijp)-dfxe
      sv(ijp) = sv(ijp)-dfye
      sw(ijp) = sw(ijp)-dfze
       
      su(ijn) = su(ijn)+dfxe
      sv(ijn) = sv(ijn)+dfye
      sw(ijn) = sw(ijn)+dfze

    end do

    ! Contribution from boundaries
    do i = 1,ninl
      iface = iInletFacesStart+i
      ijp = owner(iface)
      ijb = iInletStart + i
      call negFaceDivBoundary(arx(iface), ary(iface), arz(iface), p(ijb), su(ijp), sv(ijp), sw(ijp))
    enddo

    do i = 1,nout
      iface = iOutletFacesStart+i
      ijp = owner(iface)
      ijb = iOutletStart + i
      call negFaceDivBoundary(arx(iface), ary(iface), arz(iface), p(ijb), su(ijp), sv(ijp), sw(ijp))
    enddo

    do i = 1,nsym
      iface = iSymmetryFacesStart+i
      ijp = owner(iface)
      ijb = iSymmetryStart+i
      call negFaceDivBoundary(arx(iface), ary(iface), arz(iface), p(ijb), su(ijp), sv(ijp), sw(ijp))
    enddo   

    do i = 1,nwal
      iface = iWallFacesStart+i
      ijp = owner(iface)
      ijb = iWallStart+i
      call negFaceDivBoundary(arx(iface), ary(iface), arz(iface), p(ijb), su(ijp), sv(ijp), sw(ijp))
    enddo

    do i=1,npru
      iface = iPressOutletFacesStart + i
      ijp = owner(iface)
      ijb = iPressOutletStart + i
      call negFaceDivBoundary(arx(iface), ary(iface), arz(iface), p(ijb), su(ijp), sv(ijp), sw(ijp))
    enddo

  return
  end


!***********************************************************************
!
  function explDiv(u) result(div)
!
!***********************************************************************
!  
!  -fvc:Div(p)                                       
!  ExplDiv(u) = sum_{i=1}^{i=nf} (u)_f*sf or
!  Interpolation to cell face centers done by cds corrected scheme.
!
!***********************************************************************
!
  use geometry

  implicit none
!
!***********************************************************************
!

!...Output
    real(dp), dimension(numCells) :: div

!...Input
    real(dp), dimension(numTotal) :: u

!...Local
    integer :: i,ijp,ijn,ijb,iface
    real(dp) :: dfxe
    real(dp), dimension(3,numCells) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceDivInner(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                        u, dUdxi,dfxe)

      ! Accumulate contribution at cell center and neighbour.
      div(ijp) = div(ijp)+dfxe

      div(ijn) = div(ijn)-dfxe

    enddo

    ! Contribution from O- and C-grid cuts
    do i=1,noc
      iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
      ijp = ijl(i)
      ijn = ijr(i)
      call faceDivInner(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), &
                        u, dUdxi,dfxe)

      ! Accumulate contribution at cell center and neighbour.
      div(ijp) = div(ijp)+dfxe
       
      div(ijn) = div(ijn)-dfxe

    end do

    ! Contribution from boundaries
    do i = 1,ninl
      iface = iInletFacesStart+i
      ijp = owner(iface)
      ijb = iInletStart + i
      call FaceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), div(ijp))
    enddo

    do i = 1,nout
      iface = iOutletFacesStart+i
      ijp = owner(iface)
      ijb = iOutletStart + i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), div(ijp))
    enddo

    do i = 1,nsym
      iface = iSymmetryFacesStart+i
      ijp = owner(iface)
      ijb = iSymmetryStart+i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), div(ijp))
    enddo   

    do i = 1,nwal
      iface = iWallFacesStart+i
      ijp = owner(iface)
      ijb = iWallStart+i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), div(ijp))
    enddo

    do i=1,npru
      iface = iPressOutletFacesStart + i
      ijp = owner(iface)
      ijb = iPressOutletStart + i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), div(ijp))
    enddo

  return
  end



!
!***********************************************************************
!
  subroutine presFaceDivInner(ijp,ijn, &
                              xfc,yfc,zfc,sx,sy,sz,fif, &
                              fi,df,dfxe,dfye,dfze)
!
!***********************************************************************
! 
!     This routine calculates contribution to the explicit Divergence
!     of a scalar FI (pressure) arising from an inner cell face.
!
!***********************************************************************
!
    use types
    use parameters
    use geometry

    implicit none

    integer,    intent(in) :: ijp,ijn
    real(dp), intent(in) :: xfc,yfc,zfc
    real(dp), intent(in) :: sx,sy,sz
    real(dp), intent(in) :: fif
    real(dp), dimension(numTotal), intent(in) :: fi
    real(dp), dimension(3,numCells), intent(in) :: df


    real(dp) :: xi,yi,zi,dfxi,dfyi,dfzi
    real(dp) :: fie,dfxe,dfye,dfze
    real(dp) :: fxn,fxp
!
!***********************************************************************
!
    fxn = fif 
    fxp = 1.0d0-fxn

    xi = xc(ijp)*fxp+xc(ijn)*fxn
    yi = yc(ijp)*fxp+yc(ijn)*fxn
    zi = zc(ijp)*fxp+zc(ijn)*fxn

    dfxi = df(ijp,1)*fxp+df(ijn,1)*fxn
    dfyi = df(ijp,2)*fxp+df(ijn,2)*fxn
    dfzi = df(ijp,3)*fxp+df(ijn,3)*fxn

    ! Value of the variable at cell-face center
    fie = fi(ijp)*fxp+fi(ijn)*fxn + dfxi*(xfc-xi)+dfyi*(yfc-yi)+dfzi*(zfc-zi)

    ! (interpolated mid-face value)x(area)
    dfxe = fie*sx
    dfye = fie*sy
    dfze = fie*sz

  end subroutine


!
!***********************************************************************
!
  subroutine faceDivInner(ijp,ijn, &
                          xfc,yfc,zfc,sx,sy,sz,fif, &
                          fi,df,dfxe)
!
!***********************************************************************
! 
!     This routine calculates contribution to the explicit Divergence
!     of a scalar FI (pressure) arising from an inner cell face.
!
!***********************************************************************
!
    use types
    use parameters
    use geometry

    implicit none

    integer,    intent(in) :: ijp,ijn
    real(dp), intent(in) :: xfc,yfc,zfc
    real(dp), intent(in) :: sx,sy,sz
    real(dp), intent(in) :: fif
    real(dp), dimension(numTotal), intent(in) :: fi
    real(dp), dimension(3,numCells), intent(in) :: df


    real(dp) :: xi,yi,zi,dfxi,dfyi,dfzi
    real(dp) :: fie,dfxe
    real(dp) :: fxn,fxp
    real(dp) :: are
!
!***********************************************************************
!
    fxn = fif 
    fxp = 1.0d0-fxn

    xi = xc(ijp)*fxp+xc(ijn)*fxn
    yi = yc(ijp)*fxp+yc(ijn)*fxn
    zi = zc(ijp)*fxp+zc(ijn)*fxn

    dfxi = df(ijp,1)*fxp+df(ijn,1)*fxn
    dfyi = df(ijp,2)*fxp+df(ijn,2)*fxn
    dfzi = df(ijp,3)*fxp+df(ijn,3)*fxn

    ! Value of the variable at cell-face center
    fie = fi(ijp)*fxp+fi(ijn)*fxn + dfxi*(xfc-xi)+dfyi*(yfc-yi)+dfzi*(zfc-zi)

    ! Face area
    are = sqrt(sx**2+sy**2+sz**2)

    ! (interpolated mid-face value)x(area)
    dfxe = fie*are

  end subroutine


!***********************************************************************
!
  subroutine faceDivBoundary(sx,sy,sz,fi,dfx)
!
!***********************************************************************
!
!  This routine calculates the contribution of a boundary cell face to 
!  explicit Divergence of scalar FI.
!
!***********************************************************************
!

    implicit none

    real(dp), intent(in) :: sx,sy,sz
    real(dp), intent(in) :: fi
    real(dp), intent(inout)  :: dfx
!
!***********************************************************************
!
    dfx = dfx + fi*sqrt(sx**2+sy**2+sz**2)

  end subroutine

!***********************************************************************
!
  subroutine negFaceDivBoundary(sx,sy,sz,fi,dfx,dfy,dfz)
!
!***********************************************************************
!
!  This routine calculates the contribution of a boundary cell face to 
!  explicit Divergence of scalar FI.
!
!***********************************************************************
!

    implicit none

    real(dp), intent(in) :: sx,sy,sz
    real(dp), intent(in) :: fi
    real(dp), intent(inout)  :: dfx,dfy,dfz
!
!***********************************************************************
!
    dfx = dfx - fi*sx
    dfy = dfy - fi*sy
    dfz = dfz - fi*sz

  end subroutine



!***********************************************************************
!
  function average(u) result(aver)
!
!***********************************************************************
!                                       
!  average(u) = ( sum_{i=1}^{i=nf} (u)_f*sf) / (sum_{i=1}^{i=nf} sf)
!  Interpolation to cell face centers done by cds corrected scheme.
!
!***********************************************************************
!
  use geometry

  implicit none
!
!***********************************************************************
!

!...Output
    real(dp), dimension(numCells) :: aver

!...Input
    real(dp), dimension(numTotal) :: u

!...Local
    integer :: i,ijp,ijn,ijb,iface
    real(dp) :: dfxe
    real(dp) :: are
    real(dp), dimension(3,numCells) :: dUdxi
    real(dp), dimension(numCells) :: sumsf

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceDivInner(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), &
                        u, dUdxi,dfxe)

      ! Accumulate contribution at cell center and neighbour.
      aver(ijp) = aver(ijp)+dfxe

      aver(ijn) = aver(ijn)-dfxe

      are = sqrt( arx(i)**2 + ary(i)**2 + arz(i)**2 ) 
      sumsf(ijp) = sumsf(ijp) + are
      sumsf(ijn) = sumsf(ijn) + are

    enddo

    ! Contribution from O- and C-grid cuts
    do i=1,noc
      iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
      ijp = ijl(i)
      ijn = ijr(i)
      call faceDivInner(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), &
                        u, dUdxi,dfxe)

      ! Accumulate contribution at cell center and neighbour.
      aver(ijp) = aver(ijp)+dfxe
       
      aver(ijn) = aver(ijn)-dfxe

      are = sqrt( arx(i)**2 + ary(i)**2 + arz(i)**2 ) 
      sumsf(ijp) = sumsf(ijp) + are
      sumsf(ijn) = sumsf(ijn) + are

    end do

    ! Contribution from boundaries
    do i = 1,ninl
      iface = iInletFacesStart+i
      ijp = owner(iface)
      ijb = iInletStart + i
      call FaceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), aver(ijp))
      are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 ) 
      sumsf(ijp) = sumsf(ijp) + are
    enddo

    do i = 1,nout
      iface = iOutletFacesStart+i
      ijp = owner(iface)
      ijb = iOutletStart + i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), aver(ijp))
      are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 ) 
      sumsf(ijp) = sumsf(ijp) + are
    enddo

    do i = 1,nsym
      iface = iSymmetryFacesStart+i
      ijp = owner(iface)
      ijb = iSymmetryStart+i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), aver(ijp))
      are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 ) 
      sumsf(ijp) = sumsf(ijp) + are
    enddo   

    do i = 1,nwal
      iface = iWallFacesStart+i
      ijp = owner(iface)
      ijb = iWallStart+i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), aver(ijp))
      are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 ) 
      sumsf(ijp) = sumsf(ijp) + are
    enddo

    do i=1,npru
      iface = iPressOutletFacesStart + i
      ijp = owner(iface)
      ijb = iPressOutletStart + i
      call faceDivBoundary(arx(iface), ary(iface), arz(iface), u(ijb), aver(ijp))
      are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 ) 
      sumsf(ijp) = sumsf(ijp) + are
    enddo


    ! Divide by bounding faces total area
    do ijp = 1,numCells
      aver(ijp) = aver(ijp)/sumsf(ijp)
    enddo

  end function



!***********************************************************************
!
  function Interpolate(u) result(ui)
!
!***********************************************************************
!  
!  volScalarField -> surfaceScalarField by interpolation.
!  Interpolation to cell face centers done by cds corrected scheme.
!
!***********************************************************************
!
  use geometry

  implicit none
!
!***********************************************************************
!

!...Output
    real(dp), dimension(numFaces) :: ui

!...Input
    real(dp), dimension(numTotal) :: u

!...Local
    integer :: i,ijp,ijn,ijb,if
    real(dp), dimension(3,numCells) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceInterpolateCdsCorr( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, ui(i) )
    enddo

    ! Contribution from O- and C-grid cuts
    do i=1,noc
      if= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
      ijp = ijl(i)
      ijn = ijr(i)
      call faceInterpolateCdsCorr( ijp, ijn, xf(if), yf(if), zf(if), foc(i), u, dUdxi, ui(if) )
    end do

    ! Contribution from boundaries
    do i = 1,ninl
      if = iInletFacesStart+i
      ijb = iInletStart + i
      ui(if) = u(ijb)
    enddo

    do i = 1,nout
      if = iOutletFacesStart+i
      ijb = iOutletStart + i
      ui(if) = u(ijb)
    enddo

    do i = 1,nsym
      if = iSymmetryFacesStart+i
      ijb = iSymmetryStart+i
      ui(if) = u(ijb)
    enddo   

    do i = 1,nwal
      if = iWallFacesStart+i
      ijb = iWallStart+i
      ui(if) = u(ijb)
    enddo

    do i=1,npru
      if = iPressOutletFacesStart + i
      ijb = iPressOutletStart + i
      ui(if) = u(ijb)
    enddo

  return
  end



!***********************************************************************
!
  function surfaceSum(u) result(ssum)
!
!***********************************************************************
!  
!  volScalarField -> volScalarField by interpolation + face summation.
!  Interpolation to cell face centers done by cds corrected scheme.
!
!***********************************************************************
!
  use geometry

  implicit none
!
!***********************************************************************
!

!...Output
    real(dp), dimension(numCells) :: ssum

!...Input
    real(dp), dimension(numTotal) :: u

!...Local
    integer :: i,ijp,ijn,ijb,if
    real(dp) :: ui
    real(dp), dimension(3,numCells) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceInterpolateCdsCorr( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, ui )

      ! Accumulate contribution at cell center and neighbour.
      ssum(ijp) = ssum(ijp)+ui

      ssum(ijn) = ssum(ijn)-ui

    enddo

    ! Contribution from O- and C-grid cuts
    do i=1,noc
      if= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
      ijp = ijl(i)
      ijn = ijr(i)
      call faceInterpolateCdsCorr( ijp, ijn, xf(if), yf(if), zf(if), foc(i), u, dUdxi, ui )

      ! Accumulate contribution at cell center and neighbour.
      ssum(ijp) = ssum(ijp)+ui

      ssum(ijn) = ssum(ijn)-ui

    end do

    ! Contribution from boundaries
    do i = 1,ninl
      if = iInletFacesStart+i
      ijp = owner(if)
      ijb = iInletStart + i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo

    do i = 1,nout
      if = iOutletFacesStart+i
      ijp = owner(if)
      ijb = iOutletStart + i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo

    do i = 1,nsym
      if = iSymmetryFacesStart+i
      ijp = owner(if)
      ijb = iSymmetryStart+i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo   

    do i = 1,nwal
      if = iWallFacesStart+i
      ijp = owner(if)
      ijb = iWallStart+i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo

    do i=1,npru
      if = iPressOutletFacesStart + i
      ijp = owner(if)
      ijb = iPressOutletStart + i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo

  end function



!***********************************************************************
!
  function surfaceIntegrate(u) result(ssum)
!
!***********************************************************************
!  
!  volScalarField -> volScalarField
!  Performs a sum of face values bounding each cell and dividing by the cell volume.
!  Interpolation to cell face centers done by cds corrected scheme.
!
!***********************************************************************
!
  use geometry

  implicit none
!
!***********************************************************************
!

!...Output
    real(dp), dimension(numCells) :: ssum

!...Input
    real(dp), dimension(numTotal) :: u

!...Local
    integer :: i,ijp,ijn,ijb,if
    real(dp) :: ui
    real(dp), dimension(3,numCells) :: dUdxi

    ! Calculate cell-centered gradient
    call grad(u,dUdxi)

    ! Calculate terms integrated over surfaces

    ! Inner face
    do i=1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)
      call faceInterpolateCdsCorr( ijp, ijn, xf(i), yf(i), zf(i), facint(i), u, dUdxi, ui )

      ! Accumulate contribution at cell center and neighbour.
      ssum(ijp) = ssum(ijp)+ui

      ssum(ijn) = ssum(ijn)-ui

    enddo

    ! Contribution from O- and C-grid cuts
    do i=1,noc
      if= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
      ijp = ijl(i)
      ijn = ijr(i)
      call faceInterpolateCdsCorr( ijp, ijn, xf(if), yf(if), zf(if), foc(i), u, dUdxi, ui )

      ! Accumulate contribution at cell center and neighbour.
      ssum(ijp) = ssum(ijp)+ui

      ssum(ijn) = ssum(ijn)-ui

    end do

    ! Contribution from boundaries
    do i = 1,ninl
      if = iInletFacesStart+i
      ijp = owner(if)
      ijb = iInletStart + i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo

    do i = 1,nout
      if = iOutletFacesStart+i
      ijp = owner(if)
      ijb = iOutletStart + i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo

    do i = 1,nsym
      if = iSymmetryFacesStart+i
      ijp = owner(if)
      ijb = iSymmetryStart+i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo   

    do i = 1,nwal
      if = iWallFacesStart+i
      ijp = owner(if)
      ijb = iWallStart+i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo

    do i=1,npru
      if = iPressOutletFacesStart + i
      ijp = owner(if)
      ijb = iPressOutletStart + i
      ssum(ijp) = ssum(ijp)+u(ijb)
    enddo

    ! Divide by cell volume
    do ijp = 1,numCells
      ssum(ijp) = ssum(ijp)/vol(ijp)
    enddo

  end function


!
!***********************************************************************
!
  subroutine faceInterpolateCdsCorr(ijp,ijn, &
                                    xfc,yfc,zfc,fif, &
                                    fi,df,fie)
!
!***********************************************************************
! 
! This routine calculates face interpolant using CDS corrected.
!
!***********************************************************************
!
    use types
    use parameters
    use geometry

    implicit none

    integer,    intent(in) :: ijp,ijn
    real(dp), intent(in) :: xfc,yfc,zfc
    real(dp), intent(in) :: fif
    real(dp), dimension(numTotal), intent(in) :: fi
    real(dp), dimension(3,numCells), intent(in) :: df


    real(dp) :: xi,yi,zi,dfxi,dfyi,dfzi
    real(dp) :: fie
    real(dp) :: fxn,fxp
!
!***********************************************************************
!
    fxn = fif 
    fxp = 1.0d0-fxn

    xi = xc(ijp)*fxp+xc(ijn)*fxn
    yi = yc(ijp)*fxp+yc(ijn)*fxn
    zi = zc(ijp)*fxp+zc(ijn)*fxn

    dfxi = df(ijp,1)*fxp+df(ijn,1)*fxn
    dfyi = df(ijp,2)*fxp+df(ijn,2)*fxn
    dfzi = df(ijp,3)*fxp+df(ijn,3)*fxn

    ! Value of the variable at cell-face center
    fie = fi(ijp)*fxp+fi(ijn)*fxn + dfxi*(xfc-xi)+dfyi*(yfc-yi)+dfzi*(zfc-zi)

  end subroutine



!***********************************************************************
!
  subroutine add_random_noise_to_field(Phi,percent)
!
!***********************************************************************
!
!     Level is Max perturbation aplitude in percent (%) from given field value.
!   
!     Example usage:
!       call add_random_noise_to_field(U,10)
!
  use geometry, only: numCells,numTotal

  implicit none
!
!***********************************************************************
!
  real(dp), dimension(numTotal), intent(inout) :: phi
  integer, intent(in) :: percent
  
  integer :: inp
  real(dp) :: level
  real(dp) :: perturb

  level = dble(percent)

  do inp=1,numCells

        ! Random number based fluctuation of mean profile            
        CALL init_random_seed()
        CALL RANDOM_NUMBER(perturb)
        
        ! perturb is now between 0. and 1., we want it to be from 0 to 2*amplitude
        ! e.g. perturb = 0.9+perturb/5. when Max perturbation is +/- 10% of mean profile
        perturb = ( 1.0_dp - level/100.0_dp ) + perturb * (2*level/100.0_dp)

        Phi(INP) = perturb*Phi(inp)

  enddo

  end subroutine


! !***********************************************************************
! !
! pure function von_karman_lengthscale() result(lvk)
! !
! !***********************************************************************
! !
!   use types
!   use parameters
!   use indexes
!   use geometry
!   use variables
!   use gradients
  
!   implicit none
! !
! !***********************************************************************
! !
! !....Output
!  real(dp), dimension(numCells) :: lvk

! !....Input
!  ! (None)

! !....Locals
!  integer :: inp
!  real(dp) :: uscnd
!  real(dp) :: d2udx2,d2udy2,d2udz2
!  real(dp) :: d2vdx2,d2vdy2,d2vdz2
!  real(dp) :: d2wdx2,d2wdy2,d2wdz2
 
!   do inp=1,numCells
     
!       volr = 1./vol(inp)
      

! !.....derivatives in  x- direction: 
! !      
!       dudxe=dUdxi(1,inp+nj)*fx(inp)+dUdxi(1,inp)*(1.0d0-fx(inp))
!       dudxw=dUdxi(1,inp)*fx(inp-nj)+dUdxi(1,inp-nj)*(1.0d0-fx(inp-nj))
!       dudxn=dUdxi(1,inp+1)*fy(inp)+dUdxi(1,inp)*(1.0d0-fy(inp))
!       dudxs=dUdxi(1,inp)*fy(inp-1)+dUdxi(1,inp-1)*(1.0d0-fy(inp-1))
!       dudxt=dUdxi(1,inp+nij)*fz(inp)+dUdxi(1,inp)*(1.0d0-fz(inp))
!       dudxb=dUdxi(1,inp)*fz(inp-nij)+dUdxi(1,inp-nij)*(1.0d0-fz(inp-nij))
      
!       dvdxe=dVdxi(1,inp+nj)*fx(inp)+dVdxi(1,inp)*(1.0d0-fx(inp))
!       dvdxw=dVdxi(1,inp)*fx(inp-nj)+dVdxi(1,inp-nj)*(1.0d0-fx(inp-nj))
!       dvdxn=dVdxi(1,inp+1)*fy(inp)+dVdxi(1,inp)*(1.0d0-fy(inp))
!       dvdxs=dVdxi(1,inp)*fy(inp-1)+dVdxi(1,inp-1)*(1.0d0-fy(inp-1))
!       dvdxt=dVdxi(1,inp+nij)*fz(inp)+dVdxi(1,inp)*(1.0d0-fz(inp))
!       dvdxb=dVdxi(1,inp)*fz(inp-nij)+dVdxi(1,inp-nij)*(1.0d0-fz(inp-nij))

!       dwdxe=dWdxi(1,inp+nj)*fx(inp)+dWdxi(1,inp)*(1.0d0-fx(inp))
!       dwdxw=dWdxi(1,inp)*fx(inp-nj)+dWdxi(1,inp-nj)*(1.0d0-fx(inp-nj))
!       dwdxn=dWdxi(1,inp+1)*fy(inp)+dWdxi(1,inp)*(1.0d0-fy(inp))
!       dwdxs=dWdxi(1,inp)*fy(inp-1)+dWdxi(1,inp-1)*(1.0d0-fy(inp-1))
!       dwdxt=dWdxi(1,inp+nij)*fz(inp)+dWdxi(1,inp)*(1.0d0-fz(inp))
!       dwdxb=dWdxi(1,inp)*fz(inp-nij)+dWdxi(1,inp-nij)*(1.0d0-fz(inp-nij))
! !
! !.....derivatives in y- direction:    
!       dudye=dUdxi(2,inp+nj)*fx(inp)+dUdxi(2,inp)*(1.0d0-fx(inp))
!       dudyw=dUdxi(2,inp)*fx(inp-nj)+dUdxi(2,inp-nj)*(1.0d0-fx(inp-nj))
!       dudyn=dUdxi(2,inp+1)*fy(inp)+dUdxi(2,inp)*(1.0d0-fy(inp))
!       dudys=dUdxi(2,inp)*fy(inp-1)+dUdxi(2,inp-1)*(1.0d0-fy(inp-1))
!       dudyt=dUdxi(2,inp+nij)*fz(inp)+dUdxi(2,inp)*(1.0d0-fz(inp))
!       dudyb=dUdxi(2,inp)*fz(inp-nij)+dUdxi(2,inp-nij)*(1.0d0-fz(inp-nij))

!       dvdye=dVdxi(2,inp+nj)*fx(inp)+dVdxi(2,inp)*(1.0d0-fx(inp))
!       dvdyw=dVdxi(2,inp)*fx(inp-nj)+dVdxi(2,inp-nj)*(1.0d0-fx(inp-nj))
!       dvdyn=dVdxi(2,inp+1)*fy(inp)+dVdxi(2,inp)*(1.0d0-fy(inp))
!       dvdys=dVdxi(2,inp)*fy(inp-1)+dVdxi(2,inp-1)*(1.0d0-fy(inp-1))
!       dvdyt=dVdxi(2,inp+nij)*fz(inp)+dVdxi(2,inp)*(1.0d0-fz(inp))
!       dvdyb=dVdxi(2,inp)*fz(inp-nij)+dVdxi(2,inp-nij)*(1.0d0-fz(inp-nij))

!       dwdye=dWdxi(2,inp+nj)*fx(inp)+dWdxi(2,inp)*(1.0d0-fx(inp))
!       dwdyw=dWdxi(2,inp)*fx(inp-nj)+dWdxi(2,inp-nj)*(1.0d0-fx(inp-nj))
!       dwdyn=dWdxi(2,inp+1)*fy(inp)+dWdxi(2,inp)*(1.0d0-fy(inp))
!       dwdys=dWdxi(2,inp)*fy(inp-1)+dWdxi(2,inp-1)*(1.0d0-fy(inp-1))
!       dwdyt=dWdxi(2,inp+nij)*fz(inp)+dWdxi(2,inp)*(1.0d0-fz(inp))
!       dwdyb=dWdxi(2,inp)*fz(inp-nij)+dWdxi(2,inp-nij)*(1.0d0-fz(inp-nij))
! !
! !.....derivatives in z- direction:      
!       dudze=dUdxi(3,inp+nj)*fx(inp)+dUdxi(3,inp)*(1.0d0-fx(inp))
!       dudzw=dUdxi(3,inp)*fx(inp-nj)+dUdxi(3,inp-nj)*(1.0d0-fx(inp-nj))
!       dudzn=dUdxi(3,inp+1)*fy(inp)+dUdxi(3,inp)*(1.0d0-fy(inp))
!       dudzs=dUdxi(3,inp)*fy(inp-1)+dUdxi(3,inp-1)*(1.0d0-fy(inp-1))
!       dudzt=dUdxi(3,inp+nij)*fz(inp)+dUdxi(3,inp)*(1.0d0-fz(inp))
!       dudzb=dUdxi(3,inp)*fz(inp-nij)+dUdxi(3,inp-nij)*(1.0d0-fz(inp-nij))

!       dvdze=dVdxi(3,inp+nj)*fx(inp)+dVdxi(3,inp)*(1.0d0-fx(inp))
!       dvdzw=dVdxi(3,inp)*fx(inp-nj)+dVdxi(3,inp-nj)*(1.0d0-fx(inp-nj))
!       dvdzn=dVdxi(3,inp+1)*fy(inp)+dVdxi(3,inp)*(1.0d0-fy(inp))
!       dvdzs=dVdxi(3,inp)*fy(inp-1)+dVdxi(3,inp-1)*(1.0d0-fy(inp-1))
!       dvdzt=dVdxi(3,inp+nij)*fz(inp)+dVdxi(3,inp)*(1.0d0-fz(inp))
!       dvdzb=dVdxi(3,inp)*fz(inp-nij)+dVdxi(3,inp-nij)*(1.0d0-fz(inp-nij))

!       dwdze=dWdxi(3,inp+nj)*fx(inp)+dWdxi(3,inp)*(1.0d0-fx(inp))
!       dwdzw=dWdxi(3,inp)*fx(inp-nj)+dWdxi(3,inp-nj)*(1.0d0-fx(inp-nj))
!       dwdzn=dWdxi(3,inp+1)*fy(inp)+dWdxi(3,inp)*(1.0d0-fy(inp))
!       dwdzs=dWdxi(3,inp)*fy(inp-1)+dWdxi(3,inp-1)*(1.0d0-fy(inp-1))
!       dwdzt=dWdxi(3,inp+nij)*fz(inp)+dWdxi(3,inp)*(1.0d0-fz(inp))
!       dwdzb=dWdxi(3,inp)*fz(inp-nij)+dWdxi(3,inp-nij)*(1.0d0-fz(inp-nij))
! !
! !.....second derivatives:
!       d2udx2 = ((dudxe-dudxw)*ar1x(inp)+(dudxn-dudxs)*ar2x(inp)+ &
!                 (dudxt-dudxb)*ar3x(inp))*volr 
!       d2udy2 = ((dudye-dudyw)*ar1y(inp)+(dudyn-dudys)*ar2y(inp)+ &
!                 (dudyt-dudyb)*ar3y(inp))*volr
!       d2udz2 = ((dudze-dudzw)*ar1z(inp)+(dudzn-dudzs)*ar2z(inp)+ &
!                 (dudzt-dudzb)*ar3z(inp))*volr
! !---------------
!       d2vdx2 = ((dvdxe-dvdxw)*ar1x(inp)+(dvdxn-dvdxs)*ar2x(inp)+ &
!                 (dvdxt-dvdxb)*ar3x(inp))*volr
!       d2vdy2 = ((dvdye-dvdyw)*ar1y(inp)+(dvdyn-dvdys)*ar2y(inp)+ &
!                 (dvdyt-dvdyb)*ar3y(inp))*volr
!       d2vdz2 = ((dvdze-dvdzw)*ar1z(inp)+(dvdzn-dvdzs)*ar2z(inp)+ &
!                 (dvdzt-dvdzb)*ar3z(inp))*volr
! !---------------
!       d2wdx2 = ((dwdxe-dwdxw)*ar1x(inp)+(dwdxn-dwdxs)*ar2x(inp)+ &
!                 (dwdxt-dwdxb)*ar3x(inp))*volr
!       d2wdy2 = ((dwdye-dwdyw)*ar1y(inp)+(dwdyn-dwdys)*ar2y(inp)+ &
!                 (dwdyt-dwdyb)*ar3y(inp))*volr
!       d2wdz2 = ((dwdze-dwdzw)*ar1z(inp)+(dwdzn-dwdzs)*ar2z(inp)+ &
!                 (dwdzt-dwdzb)*ar3z(inp))*volr
! !---------------


! !...2nd velocity derivative generalized to 3d using the magnitude of
! !   velocity laplacian
!     uscnd = sqrt((d2udx2+d2udy2+d2udz2)**2+ &
!                  (d2vdx2+d2vdy2+d2vdz2)**2+ &
!                  (d2wdx2+d2wdy2+d2wdz2)**2) 

! !.....von karman length scale
!     lvk(inp) = cappa*magStrain(inp)/uscnd
                 
!   end do

! end function

    
end module fieldManipulation
