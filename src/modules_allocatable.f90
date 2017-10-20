
module types
!
! define precision for floating-point numbers
!
  ! double precision    
  integer, parameter :: dp = kind(1.0d0) 
end module types

module parameters  
  use types

  integer, parameter :: nphi=10 ! Max. no. of variables to solve
  integer, parameter :: iu=1    ! Variable identifiers
  integer, parameter :: iv=2
  integer, parameter :: iw=3
  integer, parameter :: ip=4
  integer, parameter :: ite=5
  integer, parameter :: ied=6
  integer, parameter :: ien=7
  integer, parameter :: ivis=8
  integer, parameter :: ivart=9
  integer, parameter :: icon=10

  real(dp), parameter :: one = 1.0d0
  real(dp), parameter :: zero = 0.0d0
  real(dp), parameter :: small = 1e-20
  real(dp), parameter :: great = 1e+20
  real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
  real(dp), parameter :: twothirds = 2./3._dp
  real(dp), parameter :: onethird = 1./3._dp
  
  ! Law of the wall parameters
  real(dp), parameter :: CAPPA = 0.41_dp 
  real(dp), parameter :: ELOG = 8.432_dp
  real(dp), parameter :: ctrans = 11.63


  integer :: monCell   ! Monitoring point for values, usually for log file
  integer :: pRefCell  ! Pressure reference cell
  integer :: mpoints   ! No. of monitoring points

  real(dp) :: flomas   ! mass flow at inlet
  real(dp) :: flomom   ! momentum of flow at inlet
  real(dp) :: prt1     ! Used in calcheatflux, not set at the moment TODO!
  real(dp) :: densit   ! Fluid density
  real(dp) :: viscos   ! Molecular dynamic viscosity
  real(dp) :: sormax   ! Residual toleance for SIMPLE
  real(dp) :: slarge   ! Upper limit for residuals before simulation blowup
  real(dp) :: facnap   ! Underelaxation factor for Reynolds stresses
  real(dp) :: facflx   ! Underelaxation factor for turbulent heat fluxes
  real(dp) :: uin,vin,win,tein,edin,tin,vartin,conin ! Inlet values (assumed constant accross inlet)

  real(dp) :: pranl     ! (= 0.7_dp for air, 7.0_dp for water, read it from input file.)
  logical :: lbuoy      ! Bouyancy effects - are they included in momentum and turbulence eqns. If yes we calculate heat fluxes.
  logical :: boussinesq ! Is Boussinesq hypothesis evoked yes=1/no=0, read from input file.
  real(dp) :: tref      ! Reference temperature, read from input file
  real(dp) :: beta      ! Thermal expansion coefficient, read from input file
  real(dp) :: phit      ! Parameter used for GGDH and AFM in calcheatflux subroutine, read from input file.
  real(dp) :: sksi      ! Parameter used in calcheatflux, read from input file
  real(dp) :: eta       ! Parameter used in calcheatflux, read from input file
  real(dp) :: gravx,gravy,gravz ! Components of gravity acceleration vector, read from input file
  real(dp) :: facvis            ! Under-relaxation for viscosity
  real(dp) :: erough            ! E parameter for rough wall
  real(dp) :: zzero             ! z0 - wall roughness [m] for aerodynamically rough walls


  !
  ! Timesteping control
  !
  integer :: numstep   ! Total number of timesteps
  integer :: itime     ! Current timestep index
  integer :: nzapis    ! After how many timesteps we write backup and output files
  integer :: maxit     ! Maximum number of iterations in timestep, also max. number of iterations for SIMPLE iteration
  real(dp) :: timestep ! Timestep size
  real(dp) :: time     ! Elapsed time
  real(dp) :: btime    ! Coefficient for Backward Euler time stepping algorithm, if btime = 1. => BDF2 Three time implicit (2nd order), if btime=0. Backward Euler (1st order)
  
  logical :: CoNumFix         ! Is Courant no. fixed during time-stepping
  real(dp) :: CoNumFixValue   ! Fixed value for Courant number - set in modinp for now - may be read in input
  real(dp) :: CoNum,meanCoNum ! Courant number.  


  ! Choosing discretization scheme cds, luds, smart,muscl, gamma, etc.
  ! logical :: lcds,lluds,lsmart,lavl,lmuscl,lumist...,lcds_flnt,l2nd_flnt,l2ndlim_flnt,lmuscl_flnt
  logical :: lcds = .false.
  logical :: lcdsc = .false.
  logical :: lluds = .false.
  logical :: lsmart = .false.
  logical :: lavl = .false.
  logical :: lmuscl = .false.
  logical :: lumist = .false.
  logical :: lkoren = .false.
  logical :: lcharm = .false.
  logical :: lospre = .false.
  logical :: lcds_flnt = .false.
  logical :: l2nd_flnt = .false.
  logical :: lmuscl_flnt = .false.

  logical :: flux_limiter = .false.
    
  ! Logicals, mostly read from simulation-input file:
  logical :: lturb,lread,lwrite,ltest             ! turbulent simulation, read restart file, write restart file, print residual of the linear solver,.,..      
  logical :: ltransient                           ! LTRANSIENT is TRUE for transient (non-stationary) simulations              
  logical :: levm,lasm,lles,lsgdh,lggdh,lafm      ! eddy-viscosity, algebraic stress model or LES, simple gradient or generalized gradient hypothesis, algerbaic flux model
  logical :: bdf,cn                               ! control for the time-stepping algorithm
  logical :: simple,piso,pimple                   ! control for the velocity-pressure coupling algorithm
  logical :: const_mflux                          ! control for constant flow rate 
  logical :: solveOmega, solveEpsilon, SolveTKE   ! Selfexplanatory, used in 'init'


  integer :: ncorr                   ! PISO control parameter: no. of Piso corrections.
  integer :: icorr                   ! PISO iteration no.: icorr=1..ncorr
  integer :: npcor                   ! No. of pressure-corrections; non-orthogonality correctors
  integer :: ipcorr                  ! Iteration no.: ipcorr=1..npcor
  integer :: nigrad                  ! No. of iters. for iterative cell-centered gradient calculation
  integer, parameter :: nipgrad = 2  ! No. of stages for 'pressure extrapolation at boundary + pressure gradient update' loop

  integer :: TurbModel ! Turbulence model case selector

  logical :: roughWall                                            ! Is aerodynamically rough wall assumed
  real(dp) :: magUbar, gradPcmf                                   ! Magnitude of the bulk velocity,and pressure grad that will drive the constant mass-flux flow (cmf)
  real(dp) :: sumLocalContErr, globalContErr, cumulativeContErr   ! Continuity errors

  ! Those with nphi are related to each field that we calculate U,V,W,P,T,TKE,ED...
  logical,  dimension(nphi) :: lcal      ! Logical do we calculate that particular field
  integer,  dimension(nphi) :: nsw       ! Number of allowed iterations in linear solver for each variable
  real(dp), dimension(nphi) :: sor       ! Tolerance level for linear solver iterations for each variable
  real(dp), dimension(nphi) :: resor     ! Residual
  real(dp), dimension(nphi) :: rnor      ! Residual normalisation factor
  real(dp), dimension(nphi) :: prtinv    ! Inverse Prandtl numbers, (see diffusion term discretisation)
  real(dp), dimension(nphi) :: urf       ! Underrelaxation factor
  real(dp), dimension(nphi) :: urfr      ! Recipr. value of urf: 1/urf
  real(dp), dimension(nphi) :: urfm      ! 1-urf
  real(dp), dimension(nphi) :: gds       ! Gamma blending factor [0,1] for deffered correction for convection terms: uds + gds*(uhigh-uds), 0-UDS, 1-Higher order diff.scheme

end module parameters


module hcoef
!%%%%%%%%%%%%%
  use types  
    ! Arrays used in piso algorithm 
    real(dp), dimension(:), allocatable :: h
end module hcoef



module variables
!%%%%%%%%%%%%%%
  use types

    ! These are cellwise defined variables, that is - the fields
    real(dp), dimension(:), allocatable :: u,v,w              ! Velocity components
    real(dp), dimension(:), allocatable :: flmass             ! Mass fluxes trough inner faces
    real(dp), dimension(:), allocatable :: p,pp               ! Pressure, Press. correction,  
    real(dp), dimension(:), allocatable :: te,ed              ! Turb. kin. energy, Dissipation,
    real(dp), dimension(:), allocatable :: vis                ! Effective viscosity
    real(dp), dimension(:), allocatable :: den                ! Density
    real(dp), dimension(:), allocatable :: gen                ! Turb. kin. energy generation
    real(dp), dimension(:), allocatable :: t                  ! Temperature, given in Degree Celsius.
    real(dp), dimension(:), allocatable :: vart               ! Temperature variance
    ! real(dp), dimension(:), allocatable :: edd              ! Dissipation of temperature variance
    real(dp), dimension(:), allocatable :: utt,vtt,wtt        ! Turbulent heat fluxes
    ! real(dp), dimension(:), allocatable :: ret              ! Turbulent Re number
    real(dp), dimension(:), allocatable :: con                ! Concentration
    real(dp), dimension(:), allocatable :: uu,vv,ww,uv,uw,vw  ! Reynolds stress tensor components
    real(dp), dimension(:,:), allocatable :: bij              ! Reynolds stress anisotropy tensor
    real(dp), dimension(:), allocatable :: fmi, fmo, fmoc     ! Mass fluxes trough boundary faces
    real(dp), dimension(:), allocatable :: visw,ypl           ! Effective visc. for boundary face, the y+ non-dimensional distance from wall
  
    ! values from n-1 timestep
    real(dp), dimension(:), allocatable :: uo, vo, wo
    real(dp), dimension(:), allocatable :: teo
    real(dp), dimension(:), allocatable :: edo
    real(dp), dimension(:), allocatable :: to
    real(dp), dimension(:), allocatable :: varto
    real(dp), dimension(:), allocatable :: cono

    ! values from n-2 time step
    real(dp), dimension(:), allocatable :: uoo,voo,woo
    real(dp), dimension(:), allocatable :: teoo
    real(dp), dimension(:), allocatable :: edoo
    real(dp), dimension(:), allocatable :: too
    real(dp), dimension(:), allocatable :: vartoo
    real(dp), dimension(:), allocatable :: conoo 

    ! Gradients
    real(dp),dimension(:,:), allocatable :: dUdxi,dVdxi,dWdxi
    real(dp),dimension(:,:), allocatable :: dPdxi
    real(dp),dimension(:,:), allocatable :: dTEdxi
    real(dp),dimension(:,:), allocatable :: dEDdxi
    real(dp),dimension(:,:), allocatable :: dTdxi
    real(dp),dimension(:,:), allocatable :: dCondxi
    real(dp),dimension(:,:), allocatable :: dVartdxi

    real(dp), dimension(:), allocatable :: magStrain          ! Strain magnitude
    real(dp), dimension(:), allocatable :: Vorticity          ! Vorticity magnitude

end module variables

module title_mod
!%%%%%%%%%%%%

  character(len=70) :: title
  character(len=4), dimension(10) ::  chvar = (/'  U ', '  V ', '  W ', '  P ', ' TE ', ' ED ', '  T ', ' VIS', 'VART', ' CON' /)
  character(len=7), dimension(10) ::  chvarSolver = &
  (/'U      ', 'V      ', 'W      ', 'p      ', 'k      ', 'epsilon', 'Temp   ', 'Visc   ', 'VarTemp', 'Conc   ' /)
  character(len=100):: input_file,inlet_file,grid_file,monitor_file,restart_file,out_folder_path
end module title_mod



module statistics
! Variables for collecting statistics
! these are ensamble averaged values over time steps
! used in t_rans-urans-hybrid rans/les approach
! these are saved in tecplot_stat file
  use types

    integer :: n_sample,istat,ifreq

    real(dp), dimension(:), allocatable :: u_aver,v_aver,w_aver
    real(dp), dimension(:), allocatable :: uu_aver,vv_aver,ww_aver, uv_aver,uw_aver,vw_aver
    real(dp), dimension(:), allocatable :: te_aver
    real(dp), dimension(:), allocatable :: t_aver
    real(dp), dimension(:), allocatable :: ut_aver,vt_aver,wt_aver
    real(dp), dimension(:), allocatable :: tt_aver
                                             
end module statistics