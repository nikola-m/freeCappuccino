!***********************************************************************
!
subroutine piso_using_pressure_correction
!
!***********************************************************************
!
! This implementation fo PISO algorithm follows descripton given in
! Ferziger, Peric - Computational Methods for Fluid Dynamics, 2nd ed.
! It uses PRESSURE CORRECTION as a variable, and is another way to 
! implement PISO algorithm.
!
! The algorithm is summarised in following steps, with referenced eqns. 
! from the book given in braces: 
!
!  1) Ansemble and solve momentum eq. with pressure field from previous 
!     outer iteration or timestep (Eq. 7.30)
!     Obtained velocity doesn't satisfy continuity eqn.
!     The way to fix this is to add velocity correction.
! 
!  Discussion before the next step:   
!
!     Velocity correction has two terms (Eq. 7.37)
!           ~                              
!     u'  = u'  - 1 / Apu * ( dp' /dx)   ;  v'   = ...  w' = ...
!      P     P                               P           P 
!
!     The first term is defined by Eq. 7.38, and uses velocity
!     corrections (u', v', and w') which we don't have at the momment.
!     In SIMPLE we neglect the tilde term - pretty crude.
!     In PISO we get the tilde velocities in an iterative manner:
!     In first step we neglect the tilde term, then in following steps
!     we calculate them using the freshest u', v' and w'
!     The idea is to approach to full pressure correction equation
!     Eq 7.39 in an iterative manner.
!     In first step we neglect the second term on the RHS of (7.39),
!     in following steps we don't neglect it, but calculate it using 
!     the new tildes.
!
!  2) Ansemble pressure equation (Eq. 7.39), neglecting the 2nd RHS term.
!     RHS then is divergence of 
!           m*          m*           m*
!     rho*u   ;   rho*v    ;   rho*w
!          P           P            P
!     Which means we interpolate these terms to cell face center, and 
!     dot them (perform scalar product) with cell face normal, i.e.
!     the face area vector.
!     *** Note, we now NEEED Rhie-Chow interpolation here! *** 
!     LHS matrix elements are divergence of unknown pressure correction 
!     gradient multiplied by rho/Ap, where Ap is diagonal term of momentum eq.
!     Note, divergence of gradient is Laplace operator, we can discretize
!     it in a routine that encapsulates implicit FVM calculation of  
!     Laplacian operator with rho/Ap as a coefficient.
!     
!  3) Solve the system with not so tight tolerance (rel. err. around 1e-2)
!     to get new pressure correction field in cell centers.
!     Following that calculate new pressure correction gradients.
!
!  This would be the end of one outer iteration is SIMPLE. In SIMPLE we would 
!  now correct mass flow at faces, correct velocities and pressure (using heavy 
!  under-relaxation).
!  In PISO on the other hand we continue by refining the velocity correction.
!
!  4) Find u', v', w' using Eq. (7.46), or (7.37) with only 2nd term on the RHS,
!              ~   ~   ~
!     and then u', v', w', using (7.38)
!     Add to exhisting source term the divergence in each cell of 
!         ~        ~       ~
!     rho*u', rho*v', rho*w', which is a second term in (7.39).
!     LHS of the pressure correction equation (7.39) remains the same.
!
!  5) Repeat step 3. 
!     Solve the system with not so tight tolerance (rel. err. around 1e-2)
!     to get new pressure correction field in cell centers.
!     Following that calculate new pressure correction gradients.
!
!  6) Using (7.38) and then (7.37) refine the velocity correction. In 
!     (7.38) we use vel. primes from step 4, i.e. the ones lagging in
!     PISO inner correction loop.
!
!  7) When happy with pressure correction exit the PISO loop, calculate
!     pressure correction gradients and correct velocities.
!
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use gradients
  use fieldManipulation
  use faceflux_mass
  ! use LIS_linear_solver_library


  implicit none
!
!***********************************************************************
!

  ! character(5) :: maxno
  ! character(10) :: tol
  integer :: i, k, inp, iface, ijp, ijn, istage
  real(dp) :: ppref, cap, can, fmcor


  a = 0.0_dp
  su = 0.0_dp

  ! Tentative (!) velocity gradients used for velocity interpolation: 
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

  ! > Assemble off diagonal entries of system matrix and find mass flux at faces using Rhie-Chow interpolation

  ! Internal faces:
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxmass2(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), cap, can, flmass(i))

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_value_index(i)
    a(k) = can

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_value_index(i)
    a(k) = cap

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - can

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - cap

    ! > Sources:

    su(ijp) = su(ijp) - flmass(i)
    su(ijn) = su(ijn) + flmass(i) 

  end do


  ! o- and c-grid cuts
  do i=1,noc

    iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
    ijp=ijl(i)
    ijn=ijr(i)

    call facefluxmass2(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), al(i), ar(i), fmoc(i))
    
    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - ar(i)
    
    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - al(i)

    ! > Sources:

    su(ijp) = su(ijp) - fmoc(i)
    su(ijn) = su(ijn) + fmoc(i)

  end do


  if(.not.const_mflux) call adjustMassFlow

  !+++++PISO Corrector loop+++++++++++++++++++++++++++++++++++++++++++++++++++
  do icorr=1,ncorr


  !=====Multiple pressure corrections=========================================
  do ipcorr=1,npcor

    ! Sastavi pressure_corr_laplacian, gore si vec uradio mass flux. Kad ga sastavljas rhs mora da se 
    ! obnovi sa novimkorekcijama neortogonalnosti, koje ce se nadoveati na su, sv, sw ge je massflux sacuvan.
    call pressure_corr_laplacian
   
    ! Solving pressure correction equation
    ! call dpcg(pp,ip)
    call iccg(pp,ip)
    ! call bicgstab(pp,ip) 
    ! write(maxno,'(i5)') nsw(ip)
    ! write(tol,'(es9.2)') sor(ip)
    ! write(options,'(a)') "-i cg -p ilu -ilu_fill 1 -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
    ! write(options,'(a)') "-i gmres -restart [20] -p ilut -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
    ! call solve_csr( numCells, nnz, ioffset, ja, a, su, pp )

   
    do istage=1,nipgrad

      ! Pressure corr. at boundaries (for correct calculation of pp gradient)
      call bpres(pp,istage)

      ! Calculate pressure-correction gradient and store it in pressure gradient field.
      call grad(pp,dPdxi)
  
    end do

  !======END: Multiple pressure corrections loop==============================
  enddo

    !                                                                            ~  ~   ~
    ! Initialize tilde velocities, we use temporary arrays spu, spv, sp to store u, v , w.
    if (icorr .eq. 1) then
      spu = 0.0_dp
      spv = 0.0_dp
      sp  = 0.0_dp
    endif

    !                                                                                      ~  ~   ~
    ! Get velocity corrections as in SIMPLE, we use temporary arrays spu, spv, sp to store u, v, w.
    !      
    do inp=1,numCells
        uprim(inp) = spu(inp) - dPdxi(1,inp) * vol(inp)*apu(inp)
        vprim(inp) = spv(inp) - dPdxi(2,inp) * vol(inp)*apv(inp)
        wprim(inp) = sp(inp)  - dPdxi(3,inp) * vol(inp)*apw(inp)
    enddo    


    !                                                                                            ~  ~   ~
    ! Get full velocity corrections as in Eq. 7.38 we use temporary arrays spu, spv, sp to store u, v, w.
    !    

    ! Assemble H(U) = - sum_j {a_j*U'_pj}, j - runs trough neighbour indices
    do i = 1,numInnerFaces
        ijp = owner(i)
        ijn = neighbour(i)

        k = icell_jcell_csr_value_index(i)
        spu(ijp) = spu(ijp) - h(k)*uprim(ijn) * apu(ijp)

        k = jcell_icell_csr_value_index(i)
        spu(ijn) = spu(ijn) - h(k)*uprim(ijp) * apu(ijn)

    enddo


    ! Assemble H(V) = - sum_j {a_j*V'_pj}, j - runs trough neighbour indices
    do i = 1,numInnerFaces
        ijp = owner(i)
        ijn = neighbour(i)

        k = icell_jcell_csr_value_index(i)
        spv(ijp) = spv(ijp) - h(k)*vprim(ijn) * apv(ijp)

        k = jcell_icell_csr_value_index(i)
        spv(ijn) = spv(ijn) - h(k)*vprim(ijp) * apv(ijn)

    enddo


    ! Assemble H(W) = - sum_j {a_j*W'_pj}, j - runs trough neighbour indices
    do i = 1,numInnerFaces
        ijp = owner(i)
        ijn = neighbour(i)

        k = icell_jcell_csr_value_index(i)
        sp(ijp) = sp(ijp) - h(k)*wprim(ijn) * apw(ijp)

        k = jcell_icell_csr_value_index(i)
        sp(ijn) = sp(ijn) - h(k)*wprim(ijp) * apw(ijn)

    enddo


  !Napravi novi RHS za pp jednacinu, znaci ponovo mass flow ali samo mass flux sa u, v , w i dodaj jos u~, v~, w~ koji su u spu, spv,sp
  ! I vrati gore da ponovo resava jednacinu za pp.
!  ......





  !+++++PISO Corrector loop+++++++++++++++++++++++++++++++++++++++++++++++++++
    enddo

    !.... sve dok do izvesne tacnosti nismo izracunali pp polje, polje korekcije pritiska, a 
    ! korekcije brzina su nam u uprim, vprim, wprim.
    do inp=1,numCells
        uprim(inp) = spu(inp) - dPdxi(1,inp) * vol(inp)*apu(inp)
        vprim(inp) = spv(inp) - dPdxi(2,inp) * vol(inp)*apv(inp)
        wprim(inp) = sp(inp)  - dPdxi(3,inp) * vol(inp)*apw(inp)
    enddo  

    ! Reference pressure correction - p'
    ppref = pp(pRefCell)


    !
    ! Correct mass fluxes at inner cv-faces only (only inner flux)
    !

    ! Inner faces:
    do iface=1,numInnerFaces

        ijp = owner(iface)
        ijn = neighbour(iface)

        ! (icell,jcell) matrix element:
        k = icell_jcell_csr_value_index(iface)

        flmass(iface) = flmass(iface) + a(k) * (pp(ijn)-pp(ijp))
  
    enddo

    !
    ! Correct mass fluxes at faces along O-C grid cuts.
    !
    do i=1,noc
        fmoc(i) = fmoc(i) + ar(i) * ( pp(ijr(i)) - pp(ijl(i)) )
    end do

    !
    ! Correct velocities and pressure
    !      
    do inp=1,numCells
        u(inp) = u(inp) + uprim(inp)
        v(inp) = v(inp) + vprim(inp)
        w(inp) = w(inp) + wprim(inp)
        p(inp) = p(inp) + urf(ip)*(pp(inp)-ppref)
    enddo   

    ! Explicit correction of boundary conditions 
    call correctBoundaryConditionsVelocity

    ! !.......................................................................................................!
    ! if(ipcorr.ne.npcor) then      
    ! !                                    
    ! ! The source term for the non-orthogonal corrector, also the secondary mass flux correction.
    ! !

    !   ! Clean RHS vector
    !   su = 0.0_dp

    !   do i=1,numInnerFaces                                                      
    !     ijp = owner(i)
    !     ijn = neighbour(i)

    !     call fluxmc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), fmcor)

    !     flmass(i) = flmass(i)+fmcor

    !     su(ijp) = su(ijp)-fmcor
    !     su(ijn) = su(ijn)+fmcor 

    !   enddo                                                              

    !   ! Faces along O-C grid cuts
    !   do i=1,noc
    !     iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
    !     ijp = ijl(i)
    !     ijn = ijr(i)

    !     call fluxmc(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), fmcor)

    !     fmoc(i)=fmoc(i)+fmcor

    !     su(ijp)=su(ijp)-fmcor
    !     su(ijn)=su(ijn)+fmcor

    !   end do
   
    !   write(6,'(27x,a,1pe10.3)') 'sumc  =',sum(su)

    ! !.......................................................................................................!
    ! elseif(ipcorr.eq.npcor.and.npcor.gt.1) then 
    ! !
    ! ! Non-orthogonal mass flux corrector if we reached the end of non-orthogonal correction loop.
    ! ! Why not!
    ! !

    !   ! Correct mass fluxes at inner cv-faces with second corr.                                                      
    !   do i=1,numInnerFaces                                                      
    !     ijp = owner(i)
    !     ijn = neighbour(i)

    !     call fluxmc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), fmcor)

    !     flmass(i) = flmass(i) + fmcor  

    !   enddo                                                             
                                                            
    !   ! Faces along O-C grid cuts
    !   do i=1,noc
    !     iface = iOCFacesStart+i

    !     call fluxmc(ijl(i), ijr(i), xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), fmcor)

    !     fmoc(i) = fmoc(i) + fmcor
        
    !   end do

    ! endif                                                             
    ! !.......................................................................................................!




!.....Write continuity error report:
  include 'continuityErrors.h'

end subroutine
