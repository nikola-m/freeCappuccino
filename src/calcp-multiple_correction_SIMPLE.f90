!***********************************************************************
!
subroutine calcp
!***********************************************************************
!
! Assemble and solve pressure correction equation in SIMPLE algorithm
! Correct mass fluxes, presure and velocity field
! Enables multiple pressure corrections for non-orthogonal meshes
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
  use LIS_linear_solver_library


  implicit none
!
!***********************************************************************
!

  ! character(5) :: maxno
  ! character(10) :: tol
  integer :: i, k, inp, iface, ijp, ijn, istage
  real(dp) :: sum, ppref, cap, can, fmcor


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


  ! Test continutity:
  ! if(ltest) write(6,'(20x,a,1pe10.3)') 'Initial sum  =',sum(su)



!=====Multiple pressure corrections=====================================
  do ipcorr=1,npcor

    ! Initialize pressure correction
    ! pp=0.0d0

    ! Solving pressure correction equation
    ! call dpcg(pp,ip)
    call iccg(pp,ip)
    ! call bicgstab(pp,ip) 
    ! write(maxno,'(i5)') nsw(ip)
    ! write(tol,'(es9.2)') sor(ip)
    ! write(options,'(a)') "-i cg -p ilu -ilu_fill 1 -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
    ! write(options,'(a)') "-i gmres -restart [20] -p ilut -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
    ! call solve_csr( numCells, nnz, ioffset, ja, a, su, pp )


       
    ! SECOND STEP *** CORRECTOR STAGE
   
    do istage=1,nipgrad

      ! Pressure corr. at boundaries (for correct calculation of pp gradient)
      call bpres(pp,istage)

      ! Calculate pressure-correction gradient and store it in pressure gradient field.
      call grad(pp,dPdxi)
  
    end do

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
        u(inp) = u(inp) - dPdxi(1,inp) * vol(inp)*apu(inp)
        v(inp) = v(inp) - dPdxi(2,inp) * vol(inp)*apv(inp)
        w(inp) = w(inp) - dPdxi(3,inp) * vol(inp)*apw(inp)
        p(inp) = p(inp) + urf(ip)*(pp(inp)-ppref)
    enddo   

    ! Explicit correction of boundary conditions 
    call correctBoundaryConditionsVelocity

    !.......................................................................................................!
    if(ipcorr.ne.npcor) then      
    !                                    
    ! The source term for the non-orthogonal corrector, also the secondary mass flux correction.
    !

      ! Clean RHS vector
      su = 0.0_dp

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
        iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
        ijp = ijl(i)
        ijn = ijr(i)

        call fluxmc(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), fmcor)

        fmoc(i)=fmoc(i)+fmcor

        su(ijp)=su(ijp)-fmcor
        su(ijn)=su(ijn)+fmcor

      end do
   
      ! write(6,'(27x,a,1pe10.3)') 'sumc  =',sum(su)

    !.......................................................................................................!
    elseif(ipcorr.eq.npcor.and.npcor.gt.1) then 
    !
    ! Non-orthogonal mass flux corrector if we reached the end of non-orthogonal correction loop.
    ! Why not!
    !

      ! Correct mass fluxes at inner cv-faces with second corr.                                                      
      do i=1,numInnerFaces                                                      
        ijp = owner(i)
        ijn = neighbour(i)

        call fluxmc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), fmcor)

        flmass(i) = flmass(i) + fmcor  

      enddo                                                             
                                                            
      ! Faces along O-C grid cuts
      do i=1,noc
        iface = iOCFacesStart+i

        call fluxmc(ijl(i), ijr(i), xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), fmcor)

        fmoc(i) = fmoc(i) + fmcor
        
      end do

    endif                                                             
    !.......................................................................................................!


!=END: Multiple pressure corrections loop==============================
  enddo

  ! Write continuity error report:
  include 'continuityErrors.h'

end subroutine
