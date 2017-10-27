!***********************************************************************
!
subroutine PIMPLE_multiple_correction
!
!***********************************************************************
!
! This implementation fo PISO algorithm follows descripton given in
! Ferziger, Peric - Computational Methods for Fluid Dynamics, 2nd ed.
! It uses PRESSURE instead of PRESSURE CORRECTION as a variable.
! The same approach is also used in OpenFOAM library.
! The algorithm is summarised in following steps, with referenced eqns. 
! from the book given in braces: 
!
!  1) Ansemble and solve momentum eq. with pressure field from previous 
!     outer iteration or timestep (Eq. 7.30)
!     Obtained velocity doesn't satisfy continuity eqn.
!
!  2) Find "m* with tilde" velocities at each cell center (Eq. 7.31
!     but without the last term, i.e. the pressure gradient term)
!     ~ m*    ~ m*    ~ m*
!     u       v       w
!      P       P       P
!     "These are though as velocities from which the contribution of the 
!     pressure gradient has been removed"-Ferziger,Peric
!
!     Relevant code here is get_rAU_x_UEqnH() subroutine.
!
!  3) Ansemble pressure equation (Eq. 7.35).
!     RHS is divergence of 
!         ~ m*        ~ m*         ~ m*
!     rho*u   ;   rho*v    ;   rho*w
!          P           P            P
!     Which means we interpolate these terms to cell face center, and 
!     dot them (perform scalar product) with cell face normal, i.e.
!     the face area vector.
!     Note, no need for Rhie-Chow interpolation here. 
!     LHS matrix elements are divergence of unknown pressure gradient 
!     multiplied by rho/Ap, where Ap is diagonal term from momentum eq.
!     Note, divergence of gradient is Laplace operator, we can discretize
!     it in a routine that encapsulated implicit FVM calculation of  
!     Laplacian operator with rho/Ap as a coefficient.
!
!  4) Solve the system with tight tolerance (abs error below 1e-6 I guess)
!     to get new pressure field in cell centers.
!     Following that calculate new pressure gradients in cell centers.
!
!  5) Correct "m* with tilde" velocities to get velocities that satisfy
!     continuity equation (Eq. 7.34). 
!      m    ~ m*                             ~ m*      ~ m*
!     u   = u     - 1 / Apu * ( dp /dx)   ;  v   = ...  w = ...
!      P     P                                P          P 
!
! We can continue with outer iterations in a loop until the convergence
! is achieved. This is aplicable in three ways: 1) for stationary cases 
! (under-relaxation required), 2) for nonstationary with very small timestep 
! as a noniterative time advancement method (under-relaxation not required),
! or 3) in nonstationary case with few outer iterations and some 
! under-relaxation.
! In this way they would look like SIMPLE,PISO and PIMPLE implementations
! in OpenFOAM respectively.
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use title_mod
  use gradients
  use hcoef
  use fieldmanipulation
  use faceflux_mass, only: facefluxmass_piso,fluxmc

  implicit none
!
!***********************************************************************
!
!
  integer :: i, k, inp, iface, istage
  integer :: ijp, ijn
  real(dp) :: cap, can
  ! real(dp) :: fmcor

  ! Backup coefficients array:
  h = a 

!+++++PISO Corrector loop++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO icorr=1,ncorr

    call get_rAU_x_UEqnH()

    call grad(U,dUdxi)
    call grad(V,dVdxi)
    call grad(W,dWdxi) 

    ! Initialize coefficient array and source:
    a = 0.0d0
    su = 0.0d0 

    ! > Assemble off diagonal entries of system matrix and find mass flux,
    !   accumulate diagonal entries of sysem matrix, and rhs vector stored in su array.

    ! Internal faces:
    do i = 1,numInnerFaces

      ijp = owner(i)
      ijn = neighbour(i)

      call facefluxmass_piso(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), cap, can, flmass(i))

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

      su(ijp) = su(ijp)-flmass(i)
      su(ijn) = su(ijn)+flmass(i) 

    end do


    ! O- and C-grid cuts
    do i=1,noc

      iface = iOCFacesStart + i
      ijp=ijl(i)
      ijn=ijr(i)

      call facefluxmass_piso(ijp,ijn,xf(iface),yf(iface),zf(iface),arx(iface),ary(iface),arz(iface),foc(i),al(i),ar(i),fmoc(i))


      ! > Elements on main diagonal:

      ! (icell,icell) main diagonal element
      k = diag(ijp)
      a(k) = a(k) - ar(i)

      ! (jcell,jcell) main diagonal element
      k = diag(ijn)
      a(k) = a(k) - al(i)


      ! > Sources:
      su(ijp) = su(ijp)-fmoc(i)
      su(ijn) = su(ijn)+fmoc(i)

    end do

    !// "adjusts the inlet and outlet fluxes to obey continuity, which is necessary for creating a well-posed
    !// problem where a solution for pressure exists." - Comment in OF pisoFOAM code.
    if(.not.const_mflux) call adjustMassFlow

    !!  "If you have a pressure equations with boundaries that do not fix pressure level, you have to fix a reference pressure." H.Jasak cfd-online forum
    !// In incompressible flow, only relative pressure matters.  Unless there is a pressure BC present,
    !// one cell's pressure has to be set to produce a unique pressure solution
    !     pEqn.setReference(pRefCell, pRefValue);
    !//
    a( ioffset(pRefCell):ioffset(pRefCell+1)-1 ) = 0.0_dp
    a( diag(pRefCell) ) = 1.0_dp

    ! Reference pressure
    su(pRefCell) = pp(pRefCell)


!=====Multiple pressure corrections======================================================.
    DO ipcorr=1,npcor                                                      
                                                                                         
      ! Solving pressure equation
      call iccg(pp,ip)


      ! !
      ! ! Mass flux correction and source term modification for the ipcorr-th corrector.
      ! !
      ! if(ipcorr.ne.npcor) then 

      !   ! Clean RHS vector
      !   su = 0.0d0

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
      !     iface = iOCFacesStart + i
      !     ijp = ijl(i)
      !     ijn = ijr(i)
      !     call fluxmc(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), fmcor)
      !     fmoc(i)=fmoc(i)+fmcor
      !     su(ijp)=su(ijp)-fmcor
      !     su(ijn)=su(ijn)+fmcor
      !   end do


      ! !
      ! ! We have hit the last iteration of nonorthogonality correction,
      ! ! correct mass fluxes at inner cv-faces with second correction.
      ! !
      ! elseif(ipcorr.eq.npcor.and.npcor.gt.1) then 

      !   do i=1,numInnerFaces                                                      
      !     ijp = owner(i)
      !     ijn = neighbour(i)
      !     call fluxmc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), fmcor)
      !     flmass(i) = flmass(i)+fmcor                                                                                             
      !   enddo                                                              

      !   ! Faces along O-C grid cuts
      !   do i=1,noc
      !     iface = iOCFacesStart + i
      !     ijp = ijl(i)
      !     ijn = ijr(i)
      !     call fluxmc(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), fmcor)
      !     fmoc(i)=fmoc(i)+fmcor
      !   end do  

      ! endif                                                           


!=====END:Multiple pressure corrections==================================================!
    enddo

    !
    ! Correct mass fluxes at inner cv-faces only (only inner flux)
    !

    !// On the last non-orthogonality correction, correct the flux using the most up-to-date pressure
    !// The .flux method includes contributions from all implicit terms of the pEqn (the Laplacian)
    !                    phi -= pEqn.flux();                                                                                                                                         
    ! We have hit the last iteration of nonorthogonality correction:                         

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

    ! Correct mass fluxes at faces along O-C grid cuts.
    do i=1,noc
      fmoc(i) = fmoc(i) + ar(i) * ( pp(ijr(i)) - pp(ijl(i)) )
    end do

    ! Write continuity error report:
    include 'continuityErrors.h'

    ! Explicitly relax pressure for momentum corrector
    do inp = 1,numCells
      p(inp) = p(inp) + urf(ip)*( pp(inp) - p(inp) )
    enddo
    
    do istage=1,nipgrad
      ! Pressure at boundaries.
      call bpres(p,istage)
      ! Calculate pressure gradient field.
      call grad(p,dPdxi)
    end do   

    !
    ! Correct velocities
    !      
    do inp=1,numCells
        u(inp) = u(inp) - apu(inp)*dPdxi(1,inp)*vol(inp)
        v(inp) = v(inp) - apv(inp)*dPdxi(2,inp)*vol(inp)
        w(inp) = w(inp) - apw(inp)*dPdxi(3,inp)*vol(inp)
    enddo          

    ! Explicit correction of velocity boundary conditions
    call correctBoundaryConditionsVelocity

!+++++PISO Corrector loop++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  enddo

end subroutine
