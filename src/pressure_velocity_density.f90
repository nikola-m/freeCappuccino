module compressible
!
! Contains functions for compressible flow calculations
! using peessure as referent varibale, that is using
! pressure based approach.
!  
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use gradients
  use LIS_linear_solver_library

  implicit none

public :: pressure_velocity_density_coupling

contains


!***********************************************************************
!
subroutine pressure_velocity_density_coupling
!
!***********************************************************************
!
  use faceflux_mass

  implicit none

  integer :: i, k, inp, iface, ijp, ijn, istage
  real(dp) :: sum, ppref, capd, cand, capc, canc, fmcor
  real(prec) :: C_rho
  real(prec), parameter :: RVOZD = 287.058  ! J/(kg K)

!
!***********************************************************************
!  

  ! Reset soarse matrix coefficinet values and source
  a = 0.0_dp
  su = 0.0_dp
  
  ! Volume timestepping sources
  do inp = 1, numCells

    ! We consider only BDF2 timestepping here
    C_rho = 1./(RVOZD*t(inp))
    sp(inp) = C_rho*(3*vol(inp))/(2*timestep)
    su(inp) = (3*den(inp)-4*deno(inp)+denoo(inp))*vol(inp)/(2*timestep)

  end do


 !**********************************************************************
 ! Contribution to fluxes and source as in incompressible case
 ! On the rhs we have -sum(m*_f). Mass flows at faces obtained
 ! using Rhie-Chow interpolation.
 !**********************************************************************

  ! Calculate gradients of tentative velocities, used for precise interpolation
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)


  ! > Assemble off diagonal entries of system matrix and find mass flux at faces using Rhie-Chow interpolation

  ! Internal faces:
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    call facefluxmass(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), facint(i), capd, cand, flmass(i))


    call facefluxsc( ijp, ijn, &
                     xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), &
                     flmass(i), facint(i), gam, &
                     fi, dFidxi, prtr, capc, canc, suadd )

    ! > Off-diagonal elements:

    ! (icell,jcell) matrix element:
    k = icell_jcell_csr_value_index(i)
    a(k) = cand + canc

    ! (jcell,icell) matrix element:
    k = jcell_icell_csr_value_index(i)
    a(k) = capd + capc

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - cand - canc

    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - capd - capc

    ! > Sources:

    su(ijp) = su(ijp) - flmass(i) + suadd
    su(ijn) = su(ijn) + flmass(i) - suadd 

  end do


  ! o- and c-grid cuts
  do i=1,noc

    iface= ijlFace(i) ! In the future implement Weiler-Atherton cliping algorithm to compute area vector components for non matching boundaries.
    ijp=ijl(i)
    ijn=ijr(i)

    call facefluxmass(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), foc(i), capd, cand, fmoc(i))
 

    call facefluxsc( ijp, ijn, &
                     xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), &
                     fmoc(i), foc(i), gam, &
                     srdoc(i), fi, dfidxi, prtr, capc, canc, suadd )

    ar(i) = cand + canc
    al(i) = capd + capc

    ! > Elements on main diagonal:

    ! (icell,icell) main diagonal element
    k = diag(ijp)
    a(k) = a(k) - ar(i)
    
    ! (jcell,jcell) main diagonal element
    k = diag(ijn)
    a(k) = a(k) - al(i)

    ! > Sources:

    su(ijp) = su(ijp) - fmoc(i) + suadd
    su(ijn) = su(ijn) + fmoc(i) - suadd

  end do


  if(.not.const_mflux) call adjustMassFlow


!*Multiple pressure corrections loop *******************************************************************
  do ipcorr=1,npcor

    ! Initialize pressure correction
    pp=0.0d0

    ! Solving pressure correction equation

    call bicgstab(pp,ip) 
 
    ! Using LIS solver
    ! write(maxno,'(i5)') nsw(ip)
    ! write(tol,'(es9.2)') sor(ip)
    ! write(options,'(a)') "-i gmres -restart [20] -p ilut -maxiter "//adjustl(maxno)//"-tol "//adjustl(tol)
    ! call solve_csr( numCells, nnz, ioffset, ja, a, su, pp )


       
    ! SECOND STEP *** CORRECTOR STAGE
   
    do istage=1,nipgrad

      ! Pressure corr. at boundaries (for correct calculation of pp gradient)
      call bpres(pp,istage)

      ! Calculate pressure-correction gradient and store it in pressure gradient field.
      call grad(pp,dPdxi)

    end do

    ! If simulation uses least-squares gradinets call this to get conservative pressure correction gradients.
    if ( lstsq_qr .or. lstsq_dm .or. lstsq_qr ) call grad(pp,dPdxi,'gauss_corrected')

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

        flmass(iface) = flmass(iface) + a(k) * pp(ijn)

        ! (jcell,icell) matrix element:
        k = jcell_icell_csr_value_index(i)

        flmass(iface) = flmass(iface) + a(k) * ( -pp(ijp) )
  
    enddo

    !
    ! Correct mass fluxes at faces along O-C grid cuts.
    !
    do i=1,noc
        fmoc(i) = fmoc(i) + ar(i) * pp(ijr(i)) - al(i) *  pp(ijl(i)) 
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

    endif                                                             
    !.......................................................................................................!


!*END: Multiple pressure corrections loop *******************************************************************
  enddo

  ! Write continuity error report:
  include 'continuityErrors.h'

return 
end



!***********************************************************************
!
subroutine facefluxsc(ijp, ijn, xf, yf, zf, arx, ary, arz, &
                      flmass, lambda, gam, FI, dFidxi, &
                      prtr, cap, can, suadd)
!
!***********************************************************************
!
  use types
  use parameters
  use variables, only: vis
  use interpolation

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flmass
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam 
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numCells), intent(in) :: dFidxi
  real(dp), intent(in) :: prtr
  real(dp), intent(inout) :: cap, can, suadd


! Local variables
  integer  :: nrelax
  character(len=12) :: approach
  real(dp) :: Cp,Ce
  real(dp) :: fii,fm
  real(dp) :: fcfie,fcfii,ffic

  real(dp) :: fxp,fxn

!----------------------------------------------------------------------


  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda


  ! Coeff in convection-like term
  Crho_den = 1./(RVOZD*t(inp)*den(inp))*fxp + 1./(RVOZD*t(ine)*den(ine))*fxe

  ! Upwind varijanta
  Crho_den_pj = 1./(RVOZD*t(ijn)*den(ijn))
  Crho_den_p  = 1./(RVOZD*t(ijp)*den(ijp))

  ! Convection fluxes - uds
  fm = flmass
  ce = min(fm,zero) 
  cp = max(fm,zero)

  !-------------------------------------------------------
  ! System matrix coefficients
  !-------------------------------------------------------
  cap = max(-fm,zero)*Crho_den_pj*urf(ip)
  can = max( fm,zero)*Crho_den_p *urf(ip)
  !-------------------------------------------------------


  !-------------------------------------------------------
  ! Explicit higher order convection
  !-------------------------------------------------------
  if( flmass .ge. zero ) then 
    ! Flow goes from p to pj - > p is the upwind node
    fii = face_value(ijp, ijn, xf, yf, zf, fxp, fi, dFidxi)
  else
    ! Other way, flow goes from pj, to p -> pj is the upwind node.
    fii = face_value(ijn, ijp, xf, yf, zf, fxn, fi, dFidxi)
  endif

  fcfie = fm*fii

  !-------------------------------------------------------
  ! Explicit first order convection
  !-------------------------------------------------------
  fcfii = ce*fi(ijn)+cp*fi(ijp)

  !-------------------------------------------------------
  ! Deffered correction for convection = gama_blending*(high-low)
  !-------------------------------------------------------
  ffic = gam*(fcfie-fcfii)

  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  suadd = -ffic

end subroutine




!***********************************************************************
!
subroutine facefluxsc_boundary(ijp, ijn, xf, yf, zf, arx, ary, arz, flmass, FI, dFidxi, prtr, cap, can, suadd)
!
!***********************************************************************
!
  use types
  use parameters
  use variables, only: vis

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flmass 
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numTotal), intent(in) :: dFidxi
  real(dp), intent(in) :: prtr
  real(dp), intent(inout) :: cap, can, suadd


! Local variables
  real(dp) :: Cp,Ce
  real(dp) :: fm
  real(dp) :: fxp,fxn

!----------------------------------------------------------------------


  ! Face interpolation factor
  fxn=1.0_dp
  fxp=0.0_dp

  ! Convection fluxes - uds
  fm = flmass
  ce = min(fm,zero) 
  cp = max(fm,zero)

  ! System matrix coefficients
  cap = - max(fm,zero)
  can =   min(fm,zero)

  ! if(lcds) then
  !   !---------------------------------------------
  !   ! CENTRAL DIFFERENCING SCHEME (CDS) 
  !   !---------------------------------------------
  !   ! Interpolate variable FI defined at CV centers to face using corrected CDS:
  !   !   |________Ue'___________|_______________Ucorr_____________________|
  !   fii=fi(ijp)*fxp+fi(ijn)*fxn!+dfixi*(xf-xi)+dfiyi*(yf-yi)+dfizi*(zf-zi)

  !   ! Explicit second order convection 
  !   fcfie=fm*fii
  ! else
  !   !---------------------------------------------
  !   ! Darwish-Moukalled TVD schemes for unstructured grids, IJHMT, 2003. 
  !   !---------------------------------------------
  !   ! Find r's - the gradient ratio. This is universal for all schemes.
  !   ! If flow goes from P to E
  !   r1 = (2*dFidxi(1,ijp)*xpn + 2*dFidxi(2,ijp)*ypn + 2*dFidxi(3,ijp)*zpn)/(FI(ijn)-FI(ijp)) - 1.0_dp  
  !   ! If flow goes from E to P
  !   r2 = (2*dFidxi(1,ijn)*xpn + 2*dFidxi(2,ijn)*ypn + 2*dFidxi(3,ijn)*zpn)/(FI(ijp)-FI(ijn)) - 1.0_dp 
  !   ! Find Psi for [ MUSCL ] :
  !   psiw = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
  !   psie = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
  !   ! High order flux at cell face
  !   fcfie =  ce*(fi(ijn) + fxn*psie*(fi(ijp)-fi(ijn)))+ &
  !            cp*(fi(ijp) + fxp*psiw*(fi(ijn)-fi(ijp)))
  ! endif

  ! ! Explicit first order convection
  ! fcfii = ce*fi(ijn)+cp*fi(ijp)

  ! Deffered correction for convection = gama_blending*(high-low)
  !ffic = gam*(fcfie-fcfii)

  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  ! suadd = -ffic+fdfie-fdfii da li imam explicitnu konvekciju???
  ! suadd = fdfie-fdfii 
  !-------------------------------------------------------

end subroutine


end module