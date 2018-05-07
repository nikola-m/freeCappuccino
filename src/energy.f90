module energy
!
! Implementation of sclar transport equation for mean (in the sense of Reynolds averaging) TOTAL ENTHALPY.
!
  use types
  use parameters
  use geometry
  use variables
  use scalar_fluxes, only: facefluxsc

  implicit none

  ! Constants
  real(dp), parameter :: sigt = 0.85_dp


  private 

  public :: calculate_enthalpy_field


contains


subroutine calculate_enthaply_field()
!
! Main module routine to assemble and solve energy field.
!
  use types
  use parameters
  use variables
  use gradients
  implicit none

  call calcsc(H,dHdxi,ien) ! Assemble and solve energy eq.

end subroutine



!***********************************************************************
!
subroutine calcsc(Fi,dFidxi,ifi)
!
!*********************************************************************** 
!
! Ansamble and solve transport eq.
!
!***********************************************************************
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
  real(dp) :: gam, prtr, apotime, urfrs, urfms
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: coef,dcoef,sut
  real(dp) :: fimax,fimin


! Variable specific coefficients:
  gam = gds(ifi)
  prtr = 1.0_dp/sigt

! Calculate gradient: 
  call grad(fi,dfidxi)

! Initialize coef and source arrays
  a = 0.0_dp
  su = 0.0_dp
  sp = 0.0_dp

  ! Explicit divergence of mechanical energy
  mech_en_Div_source = explDiv(flmass,K)


!
!=====================================
! VOLUME SOURCE TERMS 
!=====================================
  do inp=1,numCells


    ! Mechanical energy - Ke
    Ke = 0.5_dp*(u(inp)**2+v(inp)**2+w(inp)**2) ! 1/2 * |u|^2
    Keo = 0.5_dp*(uo(inp)**2+vo(inp)**2+wo(inp)**2) 
    Keoo = 0.5_dp*(uoo(inp)**2+voo(inp)**2+woo(inp)**2) 

    if (bdf) mech_en_ddt_source = den(inp)*vol(inp) * ( 1.5_dp*Ke - (1._dp+btime)*Keo + 0.5_dp*btime*Keoo )

    ! Explicit divergence of Ke - maybe create tmp var outside this loop and call fieldManipulation
    



    !=====================================
    ! UNSTEADY TERM
    ! Three Level Implicit Time Integration Method:
    ! in case that BTIME=0. --> Implicit Euler
    !=====================================
    if(bdf) then
      apotime=den(inp)*vol(inp)/timestep
      sut=apotime*((1._dp+btime)*to(inp)-0.5_dp*btime*too(inp))
      su(inp)=su(inp)+sut
      sp(inp)=sp(inp)+apotime*(1.0_dp+0.5_dp*btime)
    endif

    ! if(lturb.and.lbuoy) then
    !   ! When bouy activated we need the freshest utt,vtt,wtt - turbulent heat fluxes
    !   call calcheatflux 
    ! end if

  enddo



end module


+ fvc::ddt(rho, K) + fvc::div(phi, K) 