
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
  
  ! Law of the wall parameters
  real(dp), parameter :: CAPPA = 0.41_dp 
  real(dp), parameter :: ELOG = 8.432_dp
  real(dp), parameter :: ctrans = 11.63


  integer :: monCell  ! Monitoring point for values, usually for log file
  integer :: pRefCell ! Pressure reference cell
  integer :: mpoints  ! No. of monitoring points

  real(dp) :: flomas  ! mass flow at inlet
  real(dp) :: flomom  ! momentum of flow at inlet
  real(dp) :: prm1   
  real(dp) :: prt1  
  real(dp) :: densit   ! Fluid density
  real(dp) :: viscos   ! Molecular dynamic viscosity
  real(dp) :: sormax   ! Residual toleance for SIMPLE
  real(dp) :: slarge   ! Upper limit for residuals before simulation blowup
  real(dp) :: facnap   ! Underelaxation factor for Reynolds stresses
  real(dp) :: facflx   ! Underelaxation factor for turbulent heat fluxes
  real(dp) :: uin,vin,win,tein,edin,tin,vartin,conin ! Inlet values (assumed constant accross inlet)

  logical :: lbuoy
  integer :: boussinesq
  real(dp) :: tref
  real(dp) :: beta
  real(dp) :: phit
  real(dp) :: sksi
  real(dp) :: eta 
  real(dp) :: gravx,gravy,gravz
  real(dp) :: facvis

  real(dp) :: erough
  real(dp) :: zzero 

  real(dp) :: pfun
  real(dp) :: hcoef
  ! real(dp) :: t_hflx,psi,cu

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


  ! Choosing discretization scheme cds, luds, smart,muscl, gamma, more in the future...
  logical :: lcds,lluds,lsmart,lavl,lmuscl,lumist,lgamma
    
  ! Logicals, mostly read from simulation-input file:
  logical ::  lturb,lread,lwrite,ltest             ! turbulent simulation, read restart file, write restart file, print residual of the linear solver,.,..      
  logical ::  ltransient                           ! LTRANSIENT is TRUE for transient (non-stationary) simulations              
  logical ::  levm,lasm,lles,lsgdh,lggdh,lafm      ! eddy-viscosity, algebraic stress model or LES, simple gradient or generalized gradient hypothesis, algerbaic flux model
  logical ::  bdf,cn                               ! control for the time-stepping algorithm
  logical ::  simple,piso,pimple                   ! control for the velocity-pressure coupling algorithm
  logical ::  const_mflux                          ! control for constant flow rate 


  integer :: ncorr     ! PISO control parameter: no. of Piso corrections.
  integer :: icorr     ! PISO iteration no.: icorr=1..ncorr
  integer :: npcor     ! No. of pressure-corrections; non-orthogonality correctors
  integer :: ipcorr    ! Iteration no.: ipcorr=1..npcor
  integer :: nigrad    ! No. of iters. for iterative cell-centered gradient calculation
  integer, parameter :: nipgrad = 1 ! No. of stages for 'pressure extrapolation at boundary + pressure gradient update' loop

  ! Turbulence model case selector: 
  ! 1-Std k-eps,
  ! 2-k-eps with Durbin limiter,
  ! 3-RNG-k-epsilon,
  ! 4-Shih.et.al. Realizable k-epsilon,
  ! 5-k-omega Wilcox
  ! 6-k-omega Shear Stress Transport - Menter
  ! 7-EARSM - Wallin-Johansson
  ! 8-EARSM Wallin-Johansson modification by Menter
  ! 9-Scale Adaptive Simulation-SAS based on SST model
  ! 10-Seamless Alpha hybrid model based on standard k-eps model
  integer :: TurbModel ! Turbulence model case selector-see comments above (if you used 'grep -nr "TurbModel" *.f90' command)

  logical :: roughWall            ! Is aerodynamically rough wall assumed
  real(dp) :: magUbar, gradPcmf   ! Magnitude of the bulk velocity,and pressure grad that will drive the constant mass-flux flow (cmf)
  real(dp) :: sumLocalContErr, globalContErr, cumulativeContErr    ! Continuity errors

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

module OMEGA_Turb_Models
!%%%%%%%%%%%%%%
! Needed for Menter SST model and WJ-EARSM-k-omega
!%%%%%%%%%%%%%%
  use types  
    real(dp) :: sigmk1,sigmk2,sigmom1,sigmom2, & ! constants prescribed in MODINP routine
                betai1,betai2,a1,alpha1,alpha2
    real(dp) :: alpha,betta,bettast
    real(dp), dimension(:), allocatable :: domega,alphasst,bettasst, &     ! Cross diffusion coed and SST coefs
                                             qsas, &                       ! SAS model additional production term
                                             prtinv_te,prtinv_ed, &        ! 1/sigma_k; 1/sigma_epsilon
                                             cmueff                        ! size(NXYZ) ! CMUEFF the effective Cmu for EARSM
    real(dp), dimension(:,:), allocatable :: bij                           ! size(5,NXYZ) Reynolds stress anisotropy, used in EARSM
    real(dp), dimension(:), allocatable :: lvk                             ! the on Karman length scale for SAS model.
end module OMEGA_Turb_Models



module variables
!%%%%%%%%%%%%%%
  use types

    ! These are cellwise defined variables, that is - the fields
    real(dp), dimension(:), allocatable :: u,v,w  ! Velocity components
    real(dp), dimension(:), allocatable :: flmass ! Mass fluxes
    real(dp), dimension(:), allocatable :: p,pp ! Pressure, Press. correction,  
    ! real(dp), dimension(:), allocatable :: te,ed ! Turb. kin. energy, Dissipation,
    ! real(dp), dimension(:), allocatable :: t ! Temperature
    real(dp), dimension(:), allocatable :: vis ! Effective viscosity
    real(dp), dimension(:), allocatable :: den ! Density
    real(dp), dimension(:), allocatable :: gen ! turb. kin. energy generation
    ! real(dp), dimension(:), allocatable :: vart
    ! real(dp), dimension(:), allocatable :: edd
    ! real(dp), dimension(:), allocatable :: utt,vtt,wtt
    ! real(dp), dimension(:), allocatable :: ret
    ! real(dp), dimension(:), allocatable :: con
    real(dp), dimension(:), allocatable :: uu,vv,ww,uv,uw,vw
    real(dp), dimension(:), allocatable :: fmi, fmo, fmoc ! Mass fluxes trough boundary faces
    real(dp), dimension(:), allocatable :: visw,ypl       ! Effective visc. for boundary face, the y+ non-dimensional distance from wall
  
    ! values from n-1 timestep
    real(dp), dimension(:), allocatable :: uo, vo, wo
    ! real(dp), dimension(:), allocatable :: to
    ! real(dp), dimension(:), allocatable :: teo
    ! real(dp), dimension(:), allocatable :: edo
    ! real(dp), dimension(:), allocatable :: varto
    ! real(dp), dimension(:), allocatable :: cono

    ! values from n-2 time step
    real(dp), dimension(:), allocatable :: uoo,voo,woo
    ! real(dp), dimension(:), allocatable :: too
    ! real(dp), dimension(:), allocatable :: teoo
    ! real(dp), dimension(:), allocatable :: edoo
    ! real(dp), dimension(:), allocatable :: vartoo
    ! real(dp), dimension(:), allocatable :: conoo 

    real(dp),dimension(:,:), allocatable :: dUdxi,dVdxi,dWdxi
    real(dp),dimension(:,:), allocatable :: dPdxi
    ! real(dp),dimension(:,:), allocatable :: dTEdxi
    ! real(dp),dimension(:,:), allocatable :: dEDdxi
    ! real(dp),dimension(:,:), allocatable :: dTdxi
    ! real(dp),dimension(:,:), allocatable :: dCondxi
    ! real(dp),dimension(:,:), allocatable :: dVartdxi


    ! Related to seamless-alpha turbulence model:
    ! real(dp), dimension(:), allocatable :: alph,al_les,al_rans,diff    ! Related to seamless-alpha turbulence model
    ! real(dp), dimension(:), allocatable :: Timelimit                   ! Durbin Time-scale limiter
    real(dp), dimension(:), allocatable :: magStrain                   ! Strain magnitude
    real(dp), dimension(:), allocatable :: Vorticity                   ! Vorticity magnitude

end module variables

module title_mod
!%%%%%%%%%%%%

  character(len=70) :: title
  character(len=4), dimension(10) ::  chvar = (/'  U ', '  V ', '  W ', '  P ', ' TE ', ' ED ', ' IEN', ' VIS', 'VART', ' CON' /)
  character(len=7), dimension(10) ::  chvarSolver = (/'U      ', 'V      ', 'W      ', 'p      ', 'k      ', 'epsilon','Energy ',&
                                                  &   'Visc   ', 'VarTemp', 'Conc   ' /)
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
    ! real(dp), dimension(:), allocatable :: t_aver
    ! real(dp), dimension(:), allocatable :: ut_aver,vt_aver,wt_aver
    ! real(dp), dimension(:), allocatable :: tt_aver
                                             
end module statistics


