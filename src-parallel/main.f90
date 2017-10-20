!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
program cappuccino
!
!******************************************************************************
!
! Description:
!  A 3D unstructured finite volume solver for Computational Fluid Dynamics.
!  
! Usage:
! ./cappuccino <input_file> <monitor_file> <restart_file> <output_folder_path>      
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  use types
  use parameters
  use geometry
  use variables
  use title_mod
  use fieldManipulation
  use sparse_matrix, only: create_CSR_matrix_from_mesh_data,su,apu
  use utils, only: show_logo
  use gradients, only: allocate_gradients
  use temperature
  use concentration

  implicit none

  include 'mpif.h'

  integer :: iter, i, ijp, ijn, inp
  integer :: narg
  integer :: itimes, itimee
  real(dp):: magUbarStar, rUAw, gragPplus, flowDirection
  real(dp):: source
  real(dp):: suma,dt
  real :: start, finish
  character(len=9) :: timechar
!                                                                       
!******************************************************************************
!

!----------------------------------------------------------------------
! MPI start up

  call MPI_INIT(ierr)                   

  call MPI_COMM_SIZE(MPI_COMM_WORLD, & 
                               nproc,&  
                               ierr  )  

  call MPI_COMM_RANK(MPI_COMM_WORLD, & 
                               myid, &  
                               ierr  )  
  
  this = myid + 1

  if(nproc .eq. 1) then
    nproc = 0
    this = 0
  endif

  ! write(*,'(2(a,i2))') ' np = ', nproc, ' myid = ', myid
!----------------------------------------------------------------------



!  Check command line arguments
  narg=command_argument_count()
  if (narg==0.or.narg<4) write(*,'(a)') 'Usage: '&
  &'./cappuccino <input_file> <monitor_file> <restart_file> <out_folder_path>'
  call get_command_argument(1,input_file)
  call get_command_argument(2,monitor_file)
  call get_command_argument(3,restart_file)
  call get_command_argument(4,out_folder_path)

!-----------------------------------------------
!  Initialization, grid definition 
!-----------------------------------------------

  if (myid .eq. 0) then

    ! Open files
    call openfiles

    ! Print cappuccino logo to log file.
    call show_logo

    write(*,'(a)') ' '
    write(*,'(a,i2)') ' Parallel run. Number of processes, np = ', nproc
    write(*,'(a)') ' '

  endif

  ! Read input file
  call read_input_file

  ! Read and process grid files
  call mesh_geometry

  ! Allocate working arrays
  call allocate_arrays

  call allocate_gradients

  ! Create sparse matrix data structure (CSR format)
  call create_CSR_matrix_from_mesh_data

  ! Initialisation of fields
  call init

!
!===============================================
!     T i m e   l o o p : 
!===============================================

  if (myid .eq. 0) then
    write(6,'(a)') ' '
    write(6,'(a)') '  Start iteration!'
    write(6,'(a)') ' '
  endif

  itimes = itime+1
  itimee = itime+numstep

  time_loop: do itime=itimes,itimee ! time_loop 
!
!===============================================
!   Update variables : 
!===============================================

    if(bdf) then
      uoo = uo 
      voo = vo 
      woo = wo 
      teoo = teo 
      edoo = edo
      if (lcal(ien)) too = to 
      if (lcal(ivart)) vartoo = varto 
      if (lcal(icon)) conoo = cono 
    endif
    uo = u 
    vo = v 
    wo = w 
    teo = te 
    edo = ed 
    if (lcal(ien)) to = t         
    if (lcal(ivart)) varto = vart 
    if (lcal(icon)) cono = con 
!
!===============================================
!.....Set inlet boundary conditions at every timestep
!===============================================
    if(itime.eq.itimes) call bcin
!
!===============================================
!.....ITERATION CONTROL MONITOR
!===============================================
!
    ! Courant number report:
    include 'CourantNo.h'

! 
!===============================================
!.....ITERATION loop
!===============================================
!
    iteration_loop: do iter=1,maxit 

      if (myid .eq. 0) then
        call cpu_time(start)
      endif

      ! Calculate velocities.
      call calcuvw 
   
      ! Pressure-velocity coupling. Two options: SIMPLE and PISO
      if(SIMPLE)   call CALCP
      if(PISO)     call PISO_multiple_correction
      if(PIMPLE)   call PIMPLE_multiple_correction

      ! Turbulence
      if(lturb)    call correct_turbulence()

      !Scalars: Temperature , temperature variance, and concentration eqs.
      if(lcal(ien))   call calculate_temperature_field()
      ! if(lcal(ivart)) call calculate_temperature_variance_field()
      if(lcal(icon))  call calculate_concentration_field()


      if (myid .eq. 0) then
        call cpu_time(finish)
        write(timechar,'(f9.3)') finish-start
        write(6,'(3a)') '  ExecutionTime = ',trim(adjustl(timechar)),' s'
        write(6,*)
      endif


      ! ! Residual normalization, convergence check  
      ! do i=1,nphi
      !   resor(i)=resor(i)*rnor(i)
      ! end do


      ! Check residual
      source = max( resor(iu),resor(iv),resor(iw),resor(ip) ) 

      ! Now find global maximum
      call global_max( source )

      if(source.gt.slarge) then
          if ( myid .eq. 0 ) write(6,"(//,10x,a)") "*** Program terminated -  iterations diverge ***" 
          call abort_mission
      endif

      ! False time steping: skoci na next line after the end of time loop
      if(.not.ltransient) then
          if(source.lt.sormax) exit time_loop
      end if

      if(ltransient) then 

          ! Konverigao u okviru timestep-a ili potrosio sve iteracije:
          if(source.lt.sormax.or.iter.ge.maxit) then 

            if(const_mflux) then
              ! Correct driving force for a constant mass flow rate.
              include 'constant_mass_flow_forcing.f90'
            endif

            if(mod(itime,nzapis).eq.0) call writefiles
            call writehistory !<- write monitoring points
            call calc_statistics 

            ! Create and save a frame for animation
            if(mod(itime,nzapis).eq.1000) then 
               include 'create_and_save_frame.f90'
            endif

            cycle time_loop
         
          endif

      end if 

    end do iteration_loop


!===============================================
!.....Write field values after nzapis iterations 
!===============================================
    if(.not.ltransient) then
      if(mod(itime,nzapis).eq.0.and.itime.ne.numstep) then
        call writefiles
        call write_restart_files
      endif
    endif


    if(ltransient .and. myid.eq.0 ) call flush(6)

  end do time_loop

  ! False time stepping comes here after time loop with exit command

!===============================================
!.....Write field values for the next run 
!===============================================
  itime = itime-1

  call writefiles

  call write_restart_files

  ! MPI final call
  call MPI_Finalize(ierr)
      
end program