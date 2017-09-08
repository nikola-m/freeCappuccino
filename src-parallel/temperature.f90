module temperature
!
! Implementation of sclar transport equation for temperature.
!
  use types
  use parameters
  use geometry
  use variables
  use scalar_fluxes, only: facefluxsc

  implicit none

  ! Constants
  real(dp), parameter :: sigt = 0.9_dp


  private 

  public :: calculate_temperature_field


contains


subroutine calculate_temperature_field()
!
! Main module routine to assemble and solve temperature field.
!
  use types
  use parameters
  use variables
  use gradients
  implicit none

  call calcsc(T,dTdxi,ien) ! Assemble and solve temperature eq.

end subroutine



subroutine calcsc(Fi,dFidxi,ifi)
!
! Ansamble and solve transport eq. for temerature scalar.
!
  use types
  use parameters
  use geometry
  use variables
  use sparse_matrix
  use gradients
  implicit none

  integer, intent(in) :: ifi
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: dfidxi

!
! Local variables
!
  integer ::  i, k, inp, ijp, ijn, ijb, iface
  real(dp) :: gam, prtr, apotime, urfrs, urfms
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: coef,dcoef,sut
  real(dp) :: fimax,fimin


  ! Variable specific coefficients:
  gam = gds(ifi)
  prtr = 1.0d0/sigt

! Calculate gradient: 
  call grad(fi,dfidxi)

! Initialize source arrays
  su = 0.0d0
  sp = 0.0d0

!
!=====================================
! VOLUME SOURCE TERMS 
!=====================================
  do inp=1,numCells

    !
    !=====================================
    ! UNSTEADY TERM
    ! Three Level Implicit Time Integration Method:
    ! in case that BTIME=0. --> Implicit Euler
    !=====================================
    if(bdf) then
      apotime=den(inp)*vol(inp)/timestep
      sut=apotime*((1+btime)*to(inp)-0.5*btime*too(inp))
      su(inp)=su(inp)+sut
      sp(inp)=sp(inp)+apotime*(1+0.5*btime)
    endif

    if(lturb.and.lbuoy) then
      ! When bouy activated we need the freshest utt,vtt,wtt - turbulent heat fluxes
      call calcheatflux 
      call Additional_algebraic_heatflux_terms
    end if

  enddo



! Calculate terms integrated over surfaces

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
        ijp = owner(i)
        ijn = neighbour(i)

        call facefluxsc( ijp, ijn, &
                         xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), &
                         flmass(i), facint(i), gam, &
                         fi, dFidxi, prtr, cap, can, suadd, fimin, fimax )

        ! > Off-diagonal elements:

        ! (icell,jcell) matrix element:
        k = icell_jcell_csr_value_index(i)
        a(k) = can

        ! (jcell,icell) matrix element:
        k = jcell_icell_csr_value_index(i)
        a(k) = cap

        ! > Elements on main diagonal:

        ! ! (icell,icell) main diagonal element
        ! k = diag(ijp)
        ! a(k) = a(k) - can
        sp(ijp) = sp(ijp) - can

        ! ! (jcell,jcell) main diagonal element
        ! k = diag(ijn)
        ! a(k) = a(k) - cap
        sp(ijn) = sp(ijn) - cap

        ! > Sources:

        su(ijp) = su(ijp) + suadd
        su(ijn) = su(ijn) - suadd 

  enddo


  ! Contribution from o- and c-grid cuts
  do i=1,noc  
        iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
        ijp=ijl(i)
        ijn=ijr(i)

        call facefluxsc( ijp, ijn, &
                         xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                         fmoc(i), foc(i), gam, &
                         fi, dfidxi, prtr, al(i), ar(i), suadd, fimin, fimax )

        sp(ijp) = sp(ijp) - ar(i)
        sp(ijn) = sp(ijn) - al(i)

        su(ijp) = su(ijp) + suadd
        su(ijn) = su(ijn) - suadd

  end do

  ! Faces on processor boundary
  do i=1,npro
    iface = iProcFacesStart + i
    ijp = owner( iface )
    ijn = iProcStart + i

    call facefluxsc( ijp, ijn, &
                     xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                     fmpro(i), fpro(i), gam, &
                     fi, dfidxi, prtr, cap, can, suadd, fimin, fimax )

    ! > Off-diagonal elements:    
    apr(i) = can

    ! > Elements on main diagonal:
    sp(ijp) = sp(ijp) - can

    ! > Sources:
    su(ijp) = su(ijp) + suadd

  enddo

  !
  ! Boundary conditions
  !

  ! Inlet faces
  do i=1,ninl
    iface = iInletFacesStart+i
    ijp = owner(iface)
    ijb = iInletStart+i

    call facefluxsc(ijp, ijb, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), fmi(i), &
     Fi, dFidxi, prtr, cap, can, suadd)

    Sp(ijp) = Sp(ijp)-can

    Su(ijp) = Su(ijp)-can*Fi(ijb)
  end do

  ! Outlet faces
  do i=1,nout
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i

    call facefluxsc(ijp, ijb, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), fmo(i), &
     FI, dFidxi, prtr, cap, can, suadd)

    Sp(ijp) = Sp(ijp)-can

    Su(ijp) = Su(ijp)-can*Fi(ijb)
  end do

  ! Wall boundary conditions

  ! Isothermal wall boundaries

  do i=1,nwali
    iface = iWallFacesStart+i
    ijp=owner(iface)
    ijb=iWallStart+i
    dcoef = (viscos+(vis(ijp)-viscos)/sigt)/pranl ! Vrlo diskutabilno, proveriti!
    coef=dcoef*srdw(i)
    a(diag(ijp)) = a(diag(ijp)) + coef
    su(ijp) = su(ijp) + coef*t(ijb)
  end do

  ! Adiabatic wall boundaries

  do i=1,nwala
    iface = iWallFacesStart+nwali+i
    ijp=owner(iface)
    ijb=iWallStart+nwali+i
    t(ijb)=t(ijp)
  end do

  ! Modify coefficients for Crank-Nicolson
  if (cn) then

      a = 0.5_dp*a ! Doesn't affect the main diagonal because it's still zero.

        do i = 1,numInnerFaces
            ijp = owner(i)
            ijn = neighbour(i)

            k = icell_jcell_csr_value_index(i)
            su(ijp) = su(ijp) - a(k)*to(ijn)

            k = jcell_icell_csr_value_index(i)
            su(ijn) = su(ijn) - a(k)*to(ijp)
        enddo

        ! Processor boundary
        do i = 1,npro
            iface = iProcFacesStart + i 
            ijp = owner( iface )
            ijn = iProcStart + i

            ! This is relevant to previous loop over faces
            su(ijp) = su(ijp) - apr(i)*to(ijn)

            ! This is relevant to next loop over cells
            su(ijp) = su(ijp) + apr(i)*to(ijp)

        enddo  

        do ijp=1,numCells
            apotime=den(ijp)*vol(ijp)/timestep
            off_diagonal_terms = sum( a( ioffset(ijp) : ioffset(ijp+1)-1 ) ) !- a(diag(ijp))
            su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*to(ijp)
            sp(ijp) = sp(ijp)+apotime
        enddo

  endif

  ! Underrelaxation factors
  urfrs=urfr(ifi)
  urfms=urfm(ifi)

  ! Main diagonal term assembly and underrelaxation:
  do inp = 1,numCells

        ! Main diagonal term assembly:
        ! Sum all coefs in a row of a sparse matrix, but since we also included diagonal element 
        ! we substract it from the sum, to eliminate it from the sum.
        ! NOTE for parallel:
        ! Contributions to main diagonal term from neighbour cells that are in other process domain
        ! are aleady in sp at this stage.
        off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) !- a(diag(ijp)) because = 0
        a(diag(inp)) = sp(inp) - off_diagonal_terms

        ! Underelaxation:
        a(diag(inp)) = a(diag(inp))*urfrs
        su(inp) = su(inp) + urfms*a(diag(inp))*fi(inp)

  enddo

  ! Solve linear system:
  call bicgstab(fi,ifi)

!
! Update symmetry and outlet boundaries
!
  ! Symmetry faces
  do i=1,nsym
    ijp = owner(iSymmetryFacesStart+i)
    ijb = iSymmetryStart+i
    fi(ijb)=fi(ijp)
  end do
  !
 ! Outlet faces
  do i=1,nout
    ijp = owner(iOutletFacesStart+i)
    ijb = iOutletStart+i
    fi(ijb)=fi(ijp)
  end do

! These field values cannot be negative
  fi(1:numCells)=max(fi(1:numCells),small)

end subroutine calcsc

end module temperature