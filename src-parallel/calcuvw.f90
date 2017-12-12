!***********************************************************************
!
subroutine calcuvw
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use gradients, only: grad
  use faceflux_velocity, only: facefluxuvw
  use fieldManipulation, only: calcPressDiv

  implicit none
!
!***********************************************************************

!
! Local variables
!
  integer :: i, k, inp, ijp, ijn, iface
  ! integer :: istage
  real(dp) :: urfrs, urfms, apotime, heat
  real(dp) :: sut, svt, swt 
  real(dp) :: sup, svp, swp
  real(dp) :: sum_off_diagonal_terms

  logical :: ScndOrderWallBC_Model
  integer :: ijb            ! Boundary field value indexes
  real(dp) :: cp, cb        ! Temp. for eq. coeficients
  real(dp) :: cap, can      ! Temp. for eq. coeficients
  real(dp) :: vsi           ! Interpolated dynamic viscosity
  real(dp) :: cf            ! Diffusion coefficient
  real(dp) :: nxf, nyf, nzf ! Boundary face normals
  real(dp) :: Are           ! Face area
  real(dp) :: dpb           ! distance cell center to face center at boundary
  real(dp) :: vsol          ! Diffusion coefficient (II)
  real(dp) :: fdne          ! Diffusive flux auxiliary
  real(dp) :: FdUi,FdVi,FdWi! Diffusive flux auxiliary
  real(dp) :: Upb, Vpb, Wpb ! Velocity difference
  real(dp) :: Utp, Vtp, Wtp
  real(dp) :: Vnp
  real(dp) :: viss
  
  ScndOrderWallBC_Model = .false.

  ! Initialize sources
  su = 0.0_dp
  sv = 0.0_dp
  sw = 0.0_dp

  ! For u  sp => spu; for v  sp => spv; for w  sp => sp 
  spu = 0.0_dp
  spv = 0.0_dp
  sp  = 0.0_dp

  ! Velocity gradients: 
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)
  ! It can also be called with components of vector field
  ! call grad(U,V,W,dUdxi,dVdxi,dWdxi)

  ! ! Pressure gradient
  ! do istage=1,nipgrad
  !   ! Pressure at boundaries (for correct calculation of press. gradient)
  !   call bpres(p,istage)
  !   ! Calculate pressure gradient.
  !   call grad(p,dPdxi)
  ! end do

  ! Pressure divergence contribution to source
  call calcPressDiv

  ! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
  do inp=1,numCells

    !=======================================================================
    ! Pressure source terms <- instead of this we used explicit divergence above
    !=======================================================================
    ! su(inp) = -dPdxi(1,inp)*vol(inp)
    ! sv(inp) = -dPdxi(2,inp)*vol(inp)
    ! sw(inp) = -dPdxi(3,inp)*vol(inp)

    ! Constant mass flow forcing - used only on U velocity component
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(const_mflux) su(inp) = su(inp) + gradPcmf*vol(inp)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !=======================================================================
    ! Buoyancy source terms
    !=======================================================================
    if(lcal(ien).and.lbuoy) then
      !-----------------------------------------------------------------------
      !........[Boussinesq-ova aproximacija: ]
      !-----------------------------------------------------------------------
      if(boussinesq) then
        heat = beta*densit*(t(inp)-tref)*vol(inp)
      else ! if(.not.boussinesq)
        heat = (densit-den(inp))*vol(inp)
      endif
      !-----------------------------------------------------------------------
      su(inp) = su(inp) - gravx*heat
      sv(inp) = sv(inp) - gravy*heat
      sw(inp) = sw(inp) - gravz*heat
    endif

    !=======================================================================
    ! Unsteady term
    !=======================================================================
    if(bdf) then
    !-----------------------------------------------------------------------
    ! Backward differentiation formula:
    ! in case that BTIME=0. --> Implicit Euler Integration
    ! in case that BTIME=1. --> Three Level Implicit Time Integration (BDF2)
    !-----------------------------------------------------------------------
      apotime=den(inp)*vol(inp)/timestep

      sut = apotime*((1+btime)*uo(inp))
      svt = apotime*((1+btime)*vo(inp)) 
      swt = apotime*((1+btime)*wo(inp))
    
      if (btime > 0.99) then ! bdf2 scheme btime=1.
        sut = sut - apotime*(0.5*btime*uoo(inp))
        svt = svt - apotime*(0.5*btime*voo(inp))
        swt = swt - apotime*(0.5*btime*woo(inp))
      endif

      su(inp) = su(inp) + sut
      sv(inp) = sv(inp) + svt
      sw(inp) = sw(inp) + swt

      spu(inp) = spu(inp) + apotime*(1+0.5*btime)
      spv(inp) = spv(inp) + apotime*(1+0.5*btime)
      sp(inp)  = sp(inp)  + apotime*(1+0.5*btime)
    !-----------------------------------------------------------------------
    endif

  end do


  !=======================================================================
  ! Calculate Reynols stresses explicitly and additional asm terms:
  !=======================================================================
  if(lturb) then

    call calcstress

    if (lasm) call Additional_algebraic_stress_terms
    
  end if
      

  ! Calculate terms integrated over surfaces

  ! Inner faces
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxuvw( ijp, ijn, &
                      xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), &
                      flmass(i), facint(i), gds(iu), &
                      cap, can, sup, svp, swp )

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_value_index(i)
    a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_value_index(i)
    a(k) = cap

    ! > Sources: 

    su(ijp) = su(ijp) + sup
    sv(ijp) = sv(ijp) + svp
    sw(ijp) = sw(ijp) + swp

    su(ijn) = su(ijn) - sup
    sv(ijn) = sv(ijn) - svp
    sw(ijn) = sw(ijn) - swp

  end do 

  ! O- and C-grid cuts (these are not boundaries!)
  do i=1,noc
    iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
    ijp=ijl(i)
    ijn=ijr(i)

    call facefluxuvw( ijp, ijn, &
                      xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                      fmoc(i), foc(i), gds(iu), &
                      srdoc(i), al(i), ar(i), sup, svp, swp)  

    ! > Matrix coefficient contribution:
    
    spu(ijp) = spu(ijp) - ar(i)
    spv(ijp) = spv(ijp) - ar(i)
    sp(ijp)  = sp(ijp)  - ar(i)

    spu(ijn) = spu(ijn) - al(i)
    spv(ijn) = spv(ijn) - al(i)
    sp(ijn)  = sp(ijn)  - al(i)

   ! > Sources: 

    su(ijp) = su(ijp) + sup
    sv(ijp) = sv(ijp) + svp
    sw(ijp) = sw(ijp) + swp

    su(ijn) = su(ijn) - sup
    sv(ijn) = sv(ijn) - svp
    sw(ijn) = sw(ijn) - swp

  end do


  ! Faces on processor boundary
  do i=1,npro

    iface = iProcFacesStart + i
    ijp = owner( iface )
    ijn = iProcStart + i

    call facefluxuvw( ijp, ijn, &
                      xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                      fmpro(i), fpro(i), gds(iu), &
                      cap, can, sup, svp, swp )

    ! > Off-diagonal elements:    
    apr(i) = can

    ! > Matrix coefficient contribution:
    
    spu(ijp) = spu(ijp) - can
    spv(ijp) = spv(ijp) - can
    sp(ijp)  = sp(ijp)  - can

   ! > Sources: 

    su(ijp) = su(ijp) + sup
    sv(ijp) = sv(ijp) + svp
    sw(ijp) = sw(ijp) + swp

  enddo


      
  ! Implement boundary conditions


  ! Inlet (constant gradient bc)

  ! Loop Inlet faces
  do i=1,ninl
    iface = iInletFacesStart+i
    ijp = owner(iface)
    ijb = iInletStart+i

    call facefluxuvw( ijp, ijb, &
                      xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                      fmi(i), &
                      cp, cb, sup, svp, swp )

    spu(ijp) = spu(ijp) - cb
    spv(ijp) = spv(ijp) - cb
    sp(ijp)  = sp(ijp)  - cb

    su(ijp) = su(ijp) - cb*u(ijb)
    sv(ijp) = sv(ijp) - cb*v(ijb)
    sw(ijp) = sw(ijp) - cb*w(ijb)
  
  end do


  ! Outlet

  ! Loop Outlet faces
  do i=1,nout
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i

    call facefluxuvw( ijp, ijb, &
                      xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                      fmo(i), &
                      cp, cb, sup, svp, swp )  

    spu(ijp) = spu(ijp) - cb
    spv(ijp) = spv(ijp) - cb
    sp(ijp)  = sp(ijp)  - cb

    su(ijp) = su(ijp) - cb*u(ijb)
    sv(ijp) = sv(ijp) - cb*v(ijb)
    sw(ijp) = sw(ijp) - cb*w(ijb)

  end do


  ! Symmetry
  do i=1,nsym
    iface = iSymmetryFacesStart+i
    ijp = owner(iface)
    ijb = iSymmetryStart+i

    ! Diffusion coef.
    vsi = vis(ijb)
    cf = vsi*srds(i) ! cf v.s. vsol -> cf is calculated using normal distance in srds!

    ! Face area 
    are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

    ! Face normals
    nxf = arx(iface)/are
    nyf = ary(iface)/are
    nzf = arz(iface)/are

    ! Dist from cc of owner cell to cf @boundary, cannot expect dpb to be normal to boundary face in general
    dpb = sqrt( (xc(ijp)-xf(iface))**2 + (yc(ijp)-yf(iface))**2 + (zc(ijp)-zf(iface))**2 )

    ! Diffusion coef. 
    vsol = vsi*are/dpb

    ! Velocity difference vector components
    upb = u(ijp)-u(ijb)
    vpb = v(ijp)-v(ijb)
    wpb = w(ijp)-w(ijb)

    fdne = 2*cf*( upb*nxf + vpb*nyf + wpb*nzf )

    spu(ijp) = spu(ijp) + vsol
    spv(ijp) = spv(ijp) + vsol
    sp(ijp)  = sp(ijp)  + vsol

    su(ijp) = su(ijp)+vsol*u(ijp)-fdne*nxf
    sv(ijp) = sv(ijp)+vsol*v(ijp)-fdne*nyf
    sw(ijp) = sw(ijp)+vsol*w(ijp)-fdne*nzf

  end do


  do i=1,nwal
    iface = iWallFacesStart+i
    ijp = owner(iface)
    ijb = iWallStart+i

    viss = viscos ! viskoznost interolirana na boundary face
    if(lturb.and.ypl(i).gt.ctrans) viss=visw(i)
    cf=viss*srdw(i) ! cf v.s. vsol -> cf is calculated using normal distance in srdw!

    ! Face area 
    are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

    ! Face normals
    nxf = arx(iface)/are
    nyf = ary(iface)/are
    nzf = arz(iface)/are

    ! Dist from cc of owner cell to cf @boundary, cannot expect dpb to be normal to boundary face in general
    dpb = sqrt( (xc(ijp)-xf(iface))**2 + (yc(ijp)-yf(iface))**2 + (zc(ijp)-zf(iface))**2 )

    ! Diffusion coef. 
    vsol = viss*are/dpb

    ! Velocity difference vector components
    upb = u(ijp)-u(ijb)
    vpb = v(ijp)-v(ijb)
    wpb = w(ijp)-w(ijb)

    ! Velocity difference vector projected to wall face normal.
    vnp = upb*nxf+vpb*nyf+wpb*nzf

    ! Velocity difference in tangential direction.
    utp = upb-vnp*nxf
    vtp = vpb-vnp*nyf
    wtp = wpb-vnp*nzf

    if (ScndOrderWallBC_Model) then

      ! Eksplicitna difuzija
      FdUi=viss*((dUdxi(1,ijp)+dUdxi(1,ijp))*nxf+(dUdxi(2,ijp)+dVdxi(1,ijp))*nyf+(dUdxi(3,ijp)+dWdxi(1,ijp))*nzf)
      FdVi=viss*((dVdxi(1,ijp)+dUdxi(2,ijp))*nxf+(dVdxi(2,ijp)+dVdxi(2,ijp))*nyf+(dVdxi(3,ijp)+dWdxi(2,ijp))*nzf)
      FdWi=viss*((dWdxi(1,ijp)+dUdxi(3,ijp))*nxf+(dWdxi(2,ijp)+dVdxi(3,ijp))*nyf+(dWdxi(3,ijp)+dWdxi(3,ijp))*nzf)
      ! Projektujes eksplicitnu difuziju na nomalu
      FdNe = FdUi*nxf + FdVi*nyf + FdWi*nzf
      ! oduzmes od eksplicitne difuzije onu komponentu duz normale
      FdUi = FdUi-FdNe*nxf
      FdVi = FdVi-FdNe*nyf
      FdWi = FdWi-FdNe*nzf

      spu(ijp) = spu(ijp) + vsol
      spv(ijp) = spv(ijp) + vsol
      sp(ijp)  = sp(ijp)  + vsol

      Su(ijp) = Su(ijp)+Vsol*U(ijp)-(2*cf*Utp+FdUi*Are)
      Sv(ijp) = Sv(ijp)+Vsol*V(ijp)-(2*cf*Vtp+FdVi*Are)
      Sw(ijp) = Sw(ijp)+Vsol*W(ijp)-(2*cf*Wtp+FdWi*Are)

    else

      spu(ijp) = spu(ijp) + vsol
      spv(ijp) = spv(ijp) + vsol
      sp(ijp)  = sp(ijp)  + vsol

      su(ijp) = su(ijp) + vsol*u(ijp) - cf*utp
      sv(ijp) = sv(ijp) + vsol*v(ijp) - cf*vtp
      sw(ijp) = sw(ijp) + vsol*w(ijp) - cf*wtp

    endif

  enddo


  ! Modify coefficients for Crank-Nicolson
  if (cn) then
    
    a = 0.5_dp*a ! Doesn't affect the main diagonal because it's still zero.

    ! Modify coefs resulting from processor boundary faces
    apr = 0.5_dp*apr

  endif


!
!.....Assemble and solve system for U component of velocity
!

  ! Crank-Nicolson time stepping source terms
  if (cn) then

    do i = 1,numInnerFaces
        ijp = owner(i)
        ijn = neighbour(i)

        k = icell_jcell_csr_value_index(i)
        su(ijp) = su(ijp) - a(k)*uo(ijn)

        k = jcell_icell_csr_value_index(i)
        su(ijn) = su(ijn) - a(k)*uo(ijp)
    enddo

    ! Processor boundary
    do i = 1,npro
        iface = iProcFacesStart + i 
        ijp = owner( iface )
        ijn = iProcStart + i

        ! This is relevant to previous loop over faces
        su(ijp) = su(ijp) - apr(i)*uo(ijn)

        ! This is relevant to next loop over cells
        su(ijp) = su(ijp) + apr(i)*uo(ijp)

    enddo    

    do ijp=1,numCells
        apotime = den(ijp)*vol(ijp)/timestep
        sum_off_diagonal_terms = sum( a(ioffset(ijp) : ioffset(ijp+1)-1) ) - a(diag(ijp))
        su(ijp) = su(ijp) + (apotime + sum_off_diagonal_terms)*uo(ijp)
        spu(ijp) = spu(ijp) + apotime
    enddo

  endif

  urfrs=urfr(iu)
  urfms=urfm(iu)

  do inp = 1,numCells

    ! Main diagonal term assembly:
    ! Sum all coefs in a row of a sparse matrix, but since we also included diagonal element 
    ! we substract it from the sum, to eliminate it from the sum.
    ! We could also write sum( a(ioffset(inp)) : a(ioffset(inp+1)-1) ) because all diagonal terms are zero.
    ! sum_off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) - a(diag(inp)) 
    ! a(diag(inp)) = spu(inp) - sum_off_diagonal_terms

    ! NOTE for parallel:
    ! Contributions to main diagonal term from neighbour cells that are in other process domain
    ! are aleady in spu at this stage.
    a(diag(inp)) = spu(inp) 

    do k = ioffset(inp),ioffset(inp+1)-1
      if (k.eq.diag(inp)) cycle
      a(diag(inp)) = a(diag(inp)) -  a(k)
    enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = su(inp) + urfms*a(diag(inp))*u(inp)

    apu(inp) = 1./(a(diag(inp))+small)

  enddo

  ! Solve fvm equations
  call bicgstab(u,iu)

!
!.....Assemble and solve system for V component of velocity
!

  ! Crank-Nicolson time stepping source terms
  if(cn) then

    do i = 1,numInnerFaces
        ijp = owner(i)
        ijn = neighbour(i)

        k = icell_jcell_csr_value_index(i)
        sv(ijp) = sv(ijp) - a(k)*vo(ijn)

        k = jcell_icell_csr_value_index(i)
        sv(ijn) = sv(ijn) - a(k)*vo(ijp)
    enddo

    ! Processor boundary
    do i = 1,npro
        iface = iProcFacesStart + i 
        ijp = owner( iface )
        ijn = iProcStart + i

        ! This is relevant to previous loop over faces
        sv(ijp) = sv(ijp) - apr(i)*vo(ijn)

        ! This is relevant to next loop over cells
        sv(ijp) = sv(ijp) + apr(i)*vo(ijp)

    enddo

    do ijp=1,numCells
        apotime=den(ijp)*vol(ijp)/timestep
        sum_off_diagonal_terms = sum( a(ioffset(ijp) : ioffset(ijp+1)-1) ) - a(diag(ijp))
        sv(ijp) = sv(ijp) + (apotime + sum_off_diagonal_terms)*vo(ijp)
        spv(ijp) = spv(ijp)+apotime
    enddo

  endif

  urfrs=urfr(iv)
  urfms=urfm(iv)
  
  do inp=1,numCells
    a(diag(inp)) = 0.0_dp
    su(inp) = 0.0_dp
  enddo

  do inp = 1,numCells

    ! Main diagonal term assembly:
    ! sum_off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) - a(diag(inp))
    ! a(diag(inp)) = spv(inp) - sum_off_diagonal_terms

    ! NOTE for parallel:
    ! Contributions to main diagonal term from neighbour cells that are in other process domain
    ! are aleady in spv at this stage.
    a(diag(inp)) = spv(inp) 

    do k = ioffset(inp),ioffset(inp+1)-1
      if (k.eq.diag(inp)) cycle
      a(diag(inp)) = a(diag(inp)) -  a(k)
    enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = sv(inp) + urfms*a(diag(inp))*v(inp)

    apv(inp) = 1./(a(diag(inp))+small)
  enddo

  ! Solve fvm equations
  ! call jacobi(v,iv)
  call bicgstab(v,iv)

!
!.....Assemble and solve system for W component of velocity
!

  ! Crank-Nicolson time stepping source terms
  if(cn) then

    do i = 1,numInnerFaces
        ijp = owner(i)
        ijn = neighbour(i)

        k = icell_jcell_csr_value_index(i)
        sw(ijp) = sw(ijp) - a(k)*wo(ijn)

        k = jcell_icell_csr_value_index(i)
        sw(ijn) = sw(ijn) - a(k)*wo(ijp)
    enddo

    ! Processor boundary
    do i = 1,npro
        iface = iProcFacesStart + i 
        ijp = owner( iface )
        ijn = iProcStart + i

        ! This is relevant to previous loop over faces
        sw(ijp) = sw(ijp) - apr(i)*wo(ijn)

        ! This is relevant to next loop over cells
        sw(ijp) = sw(ijp) + apr(i)*wo(ijp)

    enddo

    do ijp=1,numCells
        apotime = den(ijp)*vol(ijp)/timestep
        sum_off_diagonal_terms = sum( a(ioffset(ijp) : ioffset(ijp+1)-1) ) - a(diag(ijp))
        sw(ijp) = sw(ijp) + (apotime + sum_off_diagonal_terms)*wo(ijp)
        sp(ijp) = sp(ijp) + apotime
    enddo

  endif 

  urfrs=urfr(iw)
  urfms=urfm(iw)

  do inp=1,numCells
    a(diag(inp)) = 0.0_dp
    su(inp) = 0.0_dp
  enddo
    
  do inp = 1,numCells

    ! Main diagonal term assembly:
    ! sum_off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) - a(diag(inp)) 
    ! a(diag(inp)) = sp(inp) - sum_off_diagonal_terms

    ! NOTE for parallel:
    ! Contributions to main diagonal term from neighbour cells that are in other process domain
    ! are aleady in sp at this stage.
    a(diag(inp)) = sp(inp) 

    do k = ioffset(inp),ioffset(inp+1)-1
      if (k.eq.diag(inp)) cycle
      a(diag(inp)) = a(diag(inp)) -  a(k)
    enddo

    ! Underelaxation:
    a(diag(inp)) = a(diag(inp))*urfrs
    su(inp) = sw(inp) + urfms*a(diag(inp))*w(inp)

    apw(inp) = 1./(a(diag(inp))+small)

  enddo

  ! Solve fvm equations
  call bicgstab(w,iw)

  ! MPI exchange:
  call exchange( u )
  call exchange( v )
  call exchange( w )

  ! Used in Rhie-Chow where it's interpolated to boundary, so we need it in buffer:
  call exchange( apu ) 


end subroutine