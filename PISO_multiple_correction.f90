!***********************************************************************
!
subroutine PISO_multiple_correction
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use sparse_matrix
  use variables
  use title_mod
  use gradients
  use hcoef
  use fieldmanipulation

  implicit none
!
!***********************************************************************
!
!
  integer :: i, k, inp, iface, istage
  integer :: ijp, ijn
  real(dp) :: cap, can
  real(dp) :: fmcor
  real(dp) :: sum,are

  ! Before entering the corection loop backup a_nb coefficient arrays:
  h(:) = a(:)  

!+++++PISO Corrector loop++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  DO icorr=1,ncorr

    !// From the last solution of velocity, extract the diag. term from the matrix and store the reciprocal
    !// note that the matrix coefficients are functions of U due to the non-linearity of convection.
    !            volScalarField rUA = 1.0/UEqn.A();
    !// take a Jacobi pass and update U.  See Hrv Jasak's thesis eqn. 3.137 and Henrik Rusche's thesis, eqn. 2.43
    !// UEqn.H is the right-hand side of the UEqn minus the product of (the off-diagonal terms and U).
    !// Note that since the pressure gradient is not included in the UEqn. above, 
    !// this gives us U without the pressure gradient.  Also note that UEqn.H() is a function of U.
    !   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Posle ovoga imamo novo H(u)/ap, H(v)/ap ,i H(w)/ap A.K.A. "HbyA" smesteno u U,V, i W. To je polje brzina 
    ! bez uticaja gradijenta pritiska!
    !
    call get_rAU_x_UEqnH() 
    

    a(:) = 0.0d0
    su(:) = 0.0d0 

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

      ijp=ijl(i)
      ijn=ijr(i)

      call facefluxmass_piso(ijp, ijn, xfoc(i), yfoc(i), zfoc(i), xnoc(i), ynoc(i), znoc(i), foc(i), al(i), ar(i), fmoc(i))


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


    !// adjusts the inlet and outlet fluxes to obey continuity, which is necessary for creating a well-posed
    !// problem where a solution for pressure exists.
    !     adjustPhi(phi, U, p);
    if(.not.const_mflux) call adjustMassFlow

    ! Test continutity: sum=0
    write(66,'(20x,a,1pe10.3)') ' Initial sum  =',sum(su(:))

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
                                                                                            
      ! Initialize pressure
      pp(:)=0.0d0                                                        
                                                                                         
      ! Solving pressure equation
      call iccg(pp,ip)

      !                                                                                  
      ! Mass flux correction and source term modification for the ipcorr-th corrector.
      !                
      if(ipcorr.ne.npcor) then !---------------------------------------------------------------                                             
                                                                              
        ! Clean RHS vector
        su(:) = 0.0d0

        do i=1,numInnerFaces                                                      
          ijp = owner(i)
          ijn = neighbour(i)
          call fluxmc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), fmcor)
          flmass(i) = flmass(i)+fmcor 
          su(ijp) = su(ijp)-fmcor
          su(ijn) = su(ijn)+fmcor                                                                                              
        enddo                                                              

        ! Faces along O-C grid cuts
        do i=1,noc
          ijp = ijl(i)
          ijn = ijr(i)
          call fluxmc(ijp, ijn, xfoc(i), yfoc(i), zfoc(i), xnoc(i), ynoc(i), znoc(i), foc(i), fmcor)
          fmoc(i)=fmoc(i)+fmcor
          su(ijp)=su(ijp)-fmcor
          su(ijn)=su(ijn)+fmcor
        end do

          ! Test continuity sum=0. The 'sum' should drop trough successive ipcorr corrections.
        write(66,'(20x,i1,a,/,a,1pe10.3,/,20x,a,1pe10.3)')  &
                            ipcorr,'. nonorthogonal pass:', &
                                   ' sum  =',sum(su(:)),    &
                                   '|sum| =',abs(sum(su(:)))

      endif                                                                   

                                                                                                                                                                               
                                                                                             
      !// On the last non-orthogonality correction, correct the flux using the most up-to-date pressure
      !// The .flux method includes contributions from all implicit terms of the pEqn (the Laplacian)
      !                    phi -= pEqn.flux();                                                 
                                                                                               
      ! We have hit the last iteration of nonorthogonality correction:                         
      if(icorr.eq.npcor) THEN 

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



        !
        ! Correct mass fluxes at inner cv-faces with second corr.  
        ! 

        ! Inner faces      
        do i=1,numInnerFaces                                                      
          ijp = owner(i)
          ijn = neighbour(i)
          call fluxmc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), fmcor)
          flmass(i) = flmass(i)+fmcor                                                                                             
        enddo                                                              

        ! Faces along O-C grid cuts
        do i=1,noc
          ijp = ijl(i)
          ijn = ijr(i)
          call fluxmc(ijp, ijn, xfoc(i), yfoc(i), zfoc(i), xnoc(i), ynoc(i), znoc(i), foc(i), fmcor)
          fmoc(i)=fmoc(i)+fmcor
        end do 

      endif 

    ! Write PISO continuity error report:
    include 'continuityErrors.h' 

    !=====END:Multiple pressure corrections==================================================!
    enddo



    !// Add pressure gradient to interior velocity and BC's.  Note that this pressure is not just a small
    !// correction to a previous pressure, but is the entire pressure field.  Contrast this to the use of p'
    !// in Ferziger & Peric, Eqn. 7.37.
    !// NOTE: This is whole pressure, opposite to what is done in SIMPLE: p(inp)+urf(ip)*(pp(inp)-ppref) !

    ! No under-relaxation - this is the whole pressure
    p(:) = pp(:)                                                                                                                             
      
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

!.....Explicit correction of boundary conditions 
    call correctBoundaryConditionsVelocity

    ! Write continuity error report:
    !include 'continuityErrors.h'

!+++++PISO Corrector loop++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  enddo

end subroutine
