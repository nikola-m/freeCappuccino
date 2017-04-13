module k_eqn_eddy
!
! Implementation of one-equation eddy viscosity SGS, 
! Yoshizawa et al. - Bridging between eddy-viscosity-type and second-order models using an two-scale DIA,
! 9th International symposium on turbulent shear flow, Kyoto, Japan, Vol. 3., 1993, p 23.1.1-23.1.6.
!
  use types
  use parameters
  use geometry
  use variables
  use scalar_fluxes, only: facefluxsc

  implicit none

  ! Turbulence model constants
  real(dp), parameter :: Ck = 0.094_dp
  real(dp), parameter :: Ceps = 1.048_dp

  ! Derived constants
  real(dp), parameter :: CMU25 = sqrt(sqrt(ck))
  real(dp), parameter :: CMU75 = cmu25**3


  private 

  public :: correct_turbulence_k_eqn_eddy
  public :: correct_turbulence_inlet_k_eqn_eddy


contains


subroutine correct_turbulence_k_eqn_eddy()
!
! Main module routine to solve turbulence model equations and update effective viscosity.
!
  use types
  use parameters
  use variables
  use gradients
  implicit none

  call calcsc(TE,dTEdxi,ite) ! Assemble and solve turbulence kinetic energy eq.
  call modify_mu_eff()

end subroutine



subroutine correct_turbulence_inlet_k_eqn_eddy()
!
! Update effective viscosity at inlet
!
  implicit none

  call modify_mu_eff_inlet()

end subroutine



subroutine calcsc(Fi,dFidxi,ifi)
!
!
!
  use types
  use parameters
  use geometry
  use variables
  use sparse_matrix
  use gradients
  use title_mod

  implicit none

  integer, intent(in) :: ifi
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: dfidxi

!
! Local variables
!
  integer ::  i, k, inp, ijp, ijn, ijb, iface
  real(dp) :: gam, prtr, apotime, urfrs, urfms, &
              utp, vtp, wtp, utn, vtn, wtn, &
              genp, genn, sut, &
              uttbuoy, vttbuoy, wttbuoy
  real(dp) :: cap, can, suadd
  real(dp) :: magStrainSq
  real(dp) :: off_diagonal_terms
  real(dp) :: are,nxf,nyf,nzf,vnp,xtp,ytp,ztp,ut2,tau
  real(dp) :: viss
  real(dp) :: fimax,fimin
  real(dp) :: delta,third,three_halfs

  third = 1./3.0_dp
  three_halfs = 3./2.0_dp

  ! Variable specific coefficients:
  gam = gds(ifi)

  prtr = 1.0_dp

! Calculate gradient: 
  call grad(fi,dfidxi)

! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

!
! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
!

  ! TKE volume source terms

  !=========================================================
  ! STANDARD PRODUCTION
  !=========================================================
  do inp=1,numCells
    magStrainSq = magStrain(inp)*magStrain(inp)
    gen(inp) = 2*abs(vis(inp)-viscos)*magStrainSq
  enddo

  !
  !=====================================
  ! VOLUME SOURCE TERMS 
  !=====================================
  do inp=1,numCells

    genp=max(gen(inp),zero)
    genn=min(gen(inp),zero)

    ! Add production term to the rhs:
    su(inp)=genp*vol(inp)  

    ! Add destruction term to the lhs:
    delta = vol(inp)**third
    sp(inp)=den(inp) * (Ceps*te(inp)**three_halfs/delta)* vol(inp)/(te(inp)+small)
    sp(inp)=sp(inp)-genn*vol(inp)/(te(inp)+small)


    !
    !=====================================
    ! VOLUME SOURCE TERMS: buoyancy
    !=====================================
      if(lcal(ien).and.lbuoy) then
        
        ! When bouy activated we need the freshest utt,vtt,wtt - turbulent heat fluxes
        call calcheatflux 

        if(boussinesq.eq.1) then
           uttbuoy=-gravx*den(inp)*utt(inp)*vol(inp)*beta
           vttbuoy=-gravy*den(inp)*vtt(inp)*vol(inp)*beta
           wttbuoy=-gravz*den(inp)*wtt(inp)*vol(inp)*beta
        else ! if (boussinesq.eq.0)
           uttbuoy=-gravx*den(inp)*utt(inp)*vol(inp)/(t(inp)+273.15)
           vttbuoy=-gravy*den(inp)*vtt(inp)*vol(inp)/(t(inp)+273.15)
           wttbuoy=-gravz*den(inp)*wtt(inp)*vol(inp)/(t(inp)+273.15)
        end if

        utp=max(uttbuoy,zero)
        vtp=max(vttbuoy,zero)
        wtp=max(wttbuoy,zero)
        utn=min(uttbuoy,zero)
        vtn=min(vttbuoy,zero)
        wtn=min(wttbuoy,zero)

        su(inp)=su(inp)+utp+vtp+wtp
        sp(inp)=sp(inp)-(utn+vtn+wtn)/(te(inp)+small)

      end if


      ! UNSTEADY TERM
      ! Three Level Implicit Time Integration Method:
      ! in case that BTIME=0. --> Implicit Euler
      if(bdf) then

        apotime=den(inp)*vol(inp)/timestep
        sut=apotime*(1+btime)*teo(inp)

        if (btime > 0.99) then ! bdf2 scheme btime=1.
          sut = sut - apotime*(0.5*btime*teoo(inp))
        endif

        su(inp)=su(inp)+sut
        sp(inp)=sp(inp)+apotime*(1+0.5*btime)

      endif

  ! End of TKE volume source terms
  enddo



!
! CALCULATE TERMS INTEGRATED OVER FACES
!

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxsc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), flmass(i), facint(i), gam, &
     fi, dFidxi, prtr, cap, can, suadd, fimin, fimax)

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_value_index(i)
    a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_value_index(i)
    a(k) = cap

    ! > Elements on main diagonal:

    ! ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - can

    ! ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - cap

    ! > Sources:

    su(ijp) = su(ijp) + suadd
    su(ijn) = su(ijn) - suadd 

  enddo


  ! Contribution from o- and c-grid cuts
  do i=1,noc
    iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
    ijp=ijl(i)
    ijn=ijr(i)

    call facefluxsc(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), fmoc(i), foc(i), gam, &
     fi, dfidxi, prtr, al(i), ar(i), suadd, fimin, fimax)

    sp(ijp) = sp(ijp) - ar(i)
    sp(ijn) = sp(ijn) - al(i)

    su(ijp) = su(ijp) + suadd
    su(ijn) = su(ijn) - suadd
  end do

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
    do i=1,nwal
      iface = iWallFacesStart+i
      ijp=owner(iface)
      ijb=iWallStart+i

      ! Take out standard production from wall ajdecent cells.
      su(ijp)=su(ijp)-gen(ijp)*vol(ijp) 

      viss=viscos
      if(ypl(i).gt.ctrans) viss=visw(i)

      ! Face area 
      are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

      ! Face normals
      nxf = arx(iface)/are
      nyf = ary(iface)/are
      nzf = arz(iface)/are

      ! Magnitude of a cell center velocity projected on boundary face normal
      Vnp = U(ijp)*nxf+V(ijp)*nyf+W(ijp)*nzf

      ! Tangential velocity components 
      xtp = U(ijp)-Vnp*nxf
      ytp = V(ijp)-Vnp*nyf
      ztp = W(ijp)-Vnp*nzf

      ! Its magnitude
      Vtp = xtp*xtp+ytp*ytp+ztp*ztp

      ! Tangent direction - unit vector
      xtp = xtp/vtp
      ytp = ytp/vtp
      ztp = ztp/vtp

      ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
      Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

      Tau = viss*Ut2/dnw(i)

      gen(ijp)=abs(tau)*cmu25*sqrt(te(ijp))/(dnw(i)*cappa)
      su(ijp)=su(ijp)+gen(ijp)*vol(ijp)
    end do


  ! Modify coefficients for Crank-Nicolson
  if (cn) then

      a = 0.5_dp*a ! Doesn't affect the main diagonal because it's still zero.

        do i = 1,numInnerFaces
            ijp = owner(i)
            ijn = neighbour(i)

            k = icell_jcell_csr_value_index(i)
            su(ijp) = su(ijp) - a(k)*teo(ijn)

            k = jcell_icell_csr_value_index(i)
            su(ijn) = su(ijn) - a(k)*teo(ijp)
        enddo

        do ijp=1,numCells
            apotime=den(ijp)*vol(ijp)/timestep
            off_diagonal_terms = sum( a( ioffset(ijp) : ioffset(ijp+1)-1 ) ) - a(diag(ijp))
            su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*teo(ijp)
            sp(ijp) = sp(ijp)+apotime
        enddo
        
  endif

  ! Underrelaxation factors
  urfrs=urfr(ifi)
  urfms=urfm(ifi)

  ! Main diagonal term assembly:
  do inp = 1,numCells

        ! Main diagonal term assembly:
        a(diag(inp)) = sp(inp) 
        do k = ioffset(inp),ioffset(inp+1)-1
          if (k.eq.diag(inp)) cycle
          a(diag(inp)) = a(diag(inp)) -  a(k)
        enddo

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
    iface = iSymmetryFacesStart+i
    ijp = owner(iface)
    ijb = iSymmetryStart+i
    fi(ijb)=fi(ijp)
  end do
  
  ! Outlet faces
  do i=1,nout
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i
    fi(ijb)=fi(ijp)
  end do

! Report range of scalar values and clip if negative
  fimin = minval(fi)
  fimax = maxval(fi)
  write(6,'(2x,es11.4,3a,es11.4)') fimin,' <= ',chvar(ifi),' <= ',fimax

! These field values cannot be negative
  if(fimin.lt.0.0_dp) fi = max(fi,small)

end subroutine calcsc


subroutine modify_mu_eff()
!
! Update turbulent and effective viscosity.
!
  use types
  use parameters
  use geometry
  use variables
  implicit none

  integer :: i,inp
  integer :: iface, ijp,ijb
  real(dp) :: visold
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: Ut2,Tau,Utau,viscw
  real(dp) :: third,delta

  third = 1./3.0_dp

!==============================================================================
! Loop trough cells 
!==============================================================================

  do inp=1,numCells
        ! Store old value
        visold=vis(inp)
        delta = vol(inp)**third
        ! Update effective viscosity:
        ! \mu_{eff}=\mu+\mu_t; \mu_t = C_k * k^{0.5} * Delta
        vis(inp)=viscos + den(inp)*ck*sqrt(te(inp))*delta
        ! Underelaxation
        vis(inp)=urf(ivis)*vis(inp)+(1.0_dp-urf(ivis))*visold
  enddo

!==============================================================================
! Loop trough boundary faces
!==============================================================================

  !----------------------------------------------------------------------------
  ! Wall boundaries - update Visw and Ypl
  !----------------------------------------------------------------------------
  do i=1,nwal  
    iface = iWallFacesStart+i
    ijp = owner(iface)
    ijb = iWallStart+i

    ! Face area 
    are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

    ! Face normals
    nxf = arx(iface)/are
    nyf = ary(iface)/are
    nzf = arz(iface)/are

    ! Magnitude of a cell center velocity projected on boundary face normal
    Vnp = U(ijp)*nxf+V(ijp)*nyf+W(ijp)*nzf

    ! Tangential velocity components 
    xtp = U(ijp)-Vnp*nxf
    ytp = V(ijp)-Vnp*nyf
    ztp = W(ijp)-Vnp*nzf

    ! Its magnitude
    Vtp = xtp*xtp+ytp*ytp+ztp*ztp

    ! Tangent direction
    xtp = xtp/vtp
    ytp = ytp/vtp
    ztp = ztp/vtp

    ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
    Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

    Tau = viscos*Ut2/dnw(i)
    Utau = sqrt( Tau / den(ijb) )
    ypl(i) = den(ijb)*Utau*dnw(i)/viscos

    ! ! Ima i ova varijanta u cisto turb. granicni sloj varijanti sa prvom celijom u log sloju
    ! ypl(i) = den(ijb)*cmu25*sqrt(te(ijp))*dnw(i)/viscos
    ! ! ...ovo je tehnicki receno ystar iliti y* a ne y+

    viscw = zero

    if(ypl(i) > ctrans) then
      viscw = ypl(i)*viscos*cappa/log(Elog*ypl(i))
    endif

    visw(i) = max(viscos,viscw)
    vis(ijb) = visw(i)
  enddo
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Symmetry
  !----------------------------------------------------------------------------
  do i=1,nsym
    iface = iSymmetryFacesStart+i
    ijp = owner(iface)
    ijb = iSymmetryStart+i

    Vis(ijb) = Vis(ijp)
  enddo
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Inlet
  !----------------------------------------------------------------------------
  do i=1,ninl
    iface = iInletFacesStart+i
    ijp = owner(iface)
    ijb = iInletStart+i

    Vis(ijb) = Vis(ijp)
  enddo
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Outlet
  !----------------------------------------------------------------------------
  do i=1,nout
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i

    Vis(ijb) = Vis(ijp)
  enddo
  !----------------------------------------------------------------------------

end subroutine


subroutine modify_mu_eff_inlet()
!
! Update turbulent and effective viscosity at inlet.
!
  use types
  use parameters
  use geometry,only:ninl,iInletStart,iInletFacesStart,arx,ary,arz 
  use variables

  implicit none

  integer :: i,ini,iface
  real(dp) :: delta,are

  ! Loop over inlet boundaries
  do i = 1,ninl
    iface = iInletFacesStart+i
    ini = iInletStart+i   

    are = sqrt( arx(iface)**2 + ary(iface)**2 + arz(iface)**2 )
    Delta = sqrt(are)

    ! Update effective viscosity:
    ! \mu_{eff}=\mu+\mu_t; \mu_t = C_k * k^{0.5} * Delta

    vis(ini) = viscos + den(ini)*ck*sqrt(te(ini))*delta

  enddo

  if(ninl.gt.0) then ! If case hase inlet boundaries. 
    write(*,'(4x,a)') 'mu_sgs/mu: '
    do i = 1,ninl
      ini = iInletStart+i 
      write(*,'(8x,e11.4)') vis(ini)/viscos
    enddo
  endif

end subroutine


end module k_eqn_eddy