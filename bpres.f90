!***********************************************************************
! 
subroutine bpres(p,istage)
!
!***********************************************************************
! 
! Purpose:
!   Set pressures at boundaries.
!
! Discussion:
!   istage = 1 : At boundary faces, we set the values of owner cell,
!                we need those boundary values to be able to calclulate
!                pressure or pressure correction gradients.
!   istage = 2 : We perform linear extrapolation from owner cell using 
!                previosuly calculated cell centered gradients.
!   istage = 3,etc. Same as istage=2.
!                Question is do we need more than 2 passes? 
!                Numerical experiments are needed.
!
!***********************************************************************
! 
  use types
  use parameters
  use indexes
  use variables, only: dPdxi
  use geometry

  implicit none
!
!***********************************************************************
! 
  real(dp), dimension(numTotal) :: p
  integer, intent(in) :: istage

  ! Locals:
  integer :: i, ijp, ijb
  real(dp) :: xpb, ypb, zpb

  if ( istage.eq.1 ) then

    ! Loop Boundary faces:

    ! Inlet faces
    do i=1,ninl

      ijp = owner(iInletFacesStart+i)
      ijb = iInletStart+i

      ! Takes owner cell value
      p(ijb) = p(ijp)
    end do

    ! Outlet faces
    do i=1,nout

      ijp = owner(iOutletFacesStart+i)
      ijb = iOutletStart+i

      ! Takes owner cell value
      p(ijb) = p(ijp)
    end do

    ! Symmetry faces
    do i=1,nsym

      ijp = owner(iSymmetryFacesStart+i)
      ijb = iSymmetryStart+i

      ! Takes owner cell value
      p(ijb) = p(ijp)
    end do

    ! Wall faces
    do i=1,nwal
      ijp = owner(iWallFacesStart+i)
      ijb = iWallStart+i

      ! Takes owner cell value
      p(ijb) = p(ijp)   
    end do

    ! Pressure outlet faces
    do i=1,npru
      ijp = owner(iPressOutletFacesStart+i)
      ijb = iPressOutletStart+i

      ! Takes owner cell value
      p(ijb) = p(ijp)
    end do

    ! Pressure outlet faces
    do i=1,noc
      ijp = owner(iOCFacesStart+i)
      ijb = iOCStart+i

      ! Takes owner cell value
      p(ijb) = p(ijp) 
    end do

  else ! istage==2 and higher

    ! Loop Boundary faces:

    ! Inlet faces
    do i=1,ninl

      ijp = owner(iInletFacesStart+i)
      ijb = iInletStart+i

      ! Distance vector
      xpb = xfi(i)-xc(ijp) 
      ypb = yfi(i)-yc(ijp)
      zpb = zfi(i)-zc(ijp)

      ! Linear extrapolation
      p(ijb) = p(ijp) + dPdxi(1,ijp)*xpb+dPdxi(2,ijp)*ypb+dPdxi(3,ijp)*zpb

    end do

    ! Outlet faces
    do i=1,nout

      ijp = owner(iOutletFacesStart+i)
      ijb = iOutletStart+i

      ! Distance vector
      xpb = xfo(i)-xc(ijp) 
      ypb = yfo(i)-yc(ijp)
      zpb = zfo(i)-zc(ijp)

      ! Linear extrapolation
      p(ijb) = p(ijp) + dPdxi(1,ijp)*xpb+dPdxi(2,ijp)*ypb+dPdxi(3,ijp)*zpb 

    end do

    ! Symmetry faces
    do i=1,nsym

      ijp = owner(iSymmetryFacesStart+i)
      ijb = iSymmetryStart+i

      ! Distance vector
      xpb = xfs(i)-xc(ijp) 
      ypb = yfs(i)-yc(ijp)
      zpb = zfs(i)-zc(ijp)

      ! Linear extrapolation
      p(ijb) = p(ijp) + dPdxi(1,ijp)*xpb+dPdxi(2,ijp)*ypb+dPdxi(3,ijp)*zpb
      
    end do

    ! Wall faces
    do i=1,nwal
      ijp = owner(iWallFacesStart+i)
      ijb = iWallStart+i

      ! Distance vector
      xpb = xfw(i)-xc(ijp) 
      ypb = yfw(i)-yc(ijp)
      zpb = zfw(i)-zc(ijp)

      ! Linear extrapolation
      p(ijb) = p(ijp) + dPdxi(1,ijp)*xpb+dPdxi(2,ijp)*ypb+dPdxi(3,ijp)*zpb
      
    end do

    ! Pressure outlet faces
    do i=1,npru
      ijp = owner(iPressOutletFacesStart+i)
      ijb = iPressOutletStart+i

      ! Distance vector
      xpb = xfpr(i)-xc(ijp) 
      ypb = yfpr(i)-yc(ijp)
      zpb = zfpr(i)-zc(ijp)

      ! Linear extrapolation
      p(ijb) = p(ijp) + dPdxi(1,ijp)*xpb+dPdxi(2,ijp)*ypb+dPdxi(3,ijp)*zpb
      
    end do

    ! Pressure outlet faces
    do i=1,noc
      ijp = owner(iOCFacesStart+i)
      ijb = iOCStart+i

      ! Distance vector
      xpb = xfoc(i)-xc(ijp) 
      ypb = yfoc(i)-yc(ijp)
      zpb = zfoc(i)-zc(ijp)

      ! Linear extrapolation
      p(ijb) = p(ijp) + dPdxi(1,ijp)*xpb+dPdxi(2,ijp)*ypb+dPdxi(3,ijp)*zpb
      
    end do

  endif ! istage

end subroutine
