!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!
program caffa3d
!
!***********************************************************************
!
! Description:
!  An unstructured finite volume solver.
!  
! Usage:
! ./caffa3d <input_file> <inlet_file> <monitor_file> <restart_file> <out_folder_path>      
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  use types
  use parameters
  use geometry
  use variables
  use title_mod
  use fieldManipulation
  use sparse_matrix, only: su,apu
  use k_epsilon_std, only: te,teo,teoo,ed,edo,edoo
  use temperature, only: t,to,too

  implicit none

  integer :: iter, i, ijp, ijn, inp
  real(dp):: source
  integer :: narg
  real(dp):: magUbarStar, rUAw, gragPplus, flowDirection
  real(dp):: suma,dt
  real :: start, finish
  character(len=5) :: timechar
  ! integer :: imon
  ! character(len=2) :: trpn
!                                                                       
!***********************************************************************
!

!  Check if any command line arguments are found
  narg=command_argument_count()
  if (narg==0.or.narg<5) write(*,*) 'Usage: ./caffa3d <input_file> <inlet_file> <monitor_file> <restart_file> <out_folder_path>'
  call get_command_argument(1,input_file)
  call get_command_argument(2,inlet_file)
  call get_command_argument(3,monitor_file)
  call get_command_argument(4,restart_file)
  call get_command_argument(5,out_folder_path)

!-----------------------------------------------
! Files opening
!-----------------------------------------------
  call openfiles
!
!-----------------------------------------------
!  Initialisation, grid definition 
!-----------------------------------------------
  call init

!-----------------------------------------------
!  Open files for data at monitoring points 
!-----------------------------------------------
  ! if(ltransient) then
  ! open(unit=89,file=trim(out_folder_path)//'/transient_monitoring_points')
  ! open(unit=91,file=trim(out_folder_path)//'/transient_monitoring_points_names')
  ! rewind 89
  ! rewind 91
  !  do imon=1,mpoints
  !     read(91, *) trpn
  !     open(91+imon,file=trim(out_folder_path)//"/transient_monitor_point_"//trpn,access='append')
  !     if(.not.lread) rewind(91+imon)
  !  end do
  ! end if

!
!===============================================
!     T i m e   l o o p : 
!===============================================
!
  write(6,'(a)') ' '
  write(6,'(a)') '  Start iteration!'
  write(6,'(a)') ' '

  time_loop: do itime=1,numstep ! time_loop 
!
!===============================================
!     Update variables : 
!===============================================
!

  if(bdf) then
    uoo = uo 
    voo = vo 
    woo = wo 
    ! too = to 
    teoo = teo 
    edoo = edo 
    ! vartoo = varto 
    ! conoo = cono 
  endif
    uo = u 
    vo = v 
    wo = w 
    ! to = t 
    teo = te 
    edo = ed         
    ! varto = vart 
    ! cono = con 
!
!===============================================
!.....Set inlet boundary conditions at every timestep
!===============================================
      if(itime.eq.1) call bcin ! here called only once at the beginning
!
!===============================================
!.....ITERATION CONTROL MONITOR
!===============================================
!
      ! Courant number report:
      include 'CourantNo.h'

      ! Old style:
      ! Header written to monitor at start of SIMPLE iterations 
      !include 'simpleMonitorHeader.h'
      !time = time + timestep
      !write(6,*)
      !write(6,'(a,i0,a,f12.6)') ' Time step no. : ',ITIME,' Time = ',TIME
      !write(6,*)

! 
!===============================================
!.....ITERATION loop
!===============================================
!
      iteration_loop: do iter=1,maxit 

      call cpu_time(start)

!.....Calculate velocities.
      call calcuvw

!.....Update OUTLET BC.
      ! // IF(.not.const_mflux)    CALL OUTBC  

!.....Pressure-velocity coupling. Two options: SIMPLE and PISO
      IF(SIMPLE)   CALL CALCP
      IF(PISO)     CALL PISO_multiple_correction
      IF(PIMPLE)   CALL PIMPLE_multiple_correction

!.....Turbulence
      call correct_turbulence()

! !.....Turbulence equations
!       IF(LCAL(ITE))   CALL CALCSCM(TE,dTEdxi  ,ITE) 
!       IF(LCAL(IED))   CALL CALCSCM(ED,dEDdxi  ,IED)
! !.....Update turbulent eddy viscosity
!       IF( LCAL(IVIS) .AND. .NOT. LLES)  CALL MODVIS
!       ! The SGS viscosity of les model
!       IF(LCAL(IVIS) .AND. LLES)         CALL CALC_VIS_LES

! !.....Temperature , temperature variance, and concentration eqs.
!       IF(LCAL(IEN))   CALL CALCSCT(T,   dTdxi     ,IEN)
!       IF(LCAL(IVART)) CALL CALCSCT(VART,dVartdxi  ,IVART)
!       IF(LCAL(ICON))  CALL CALCSCT(CON, dCondxi   ,ICON)


      call cpu_time(finish)
      write(6,'(a,f6.3,a)') '  ExecutionTime = ',finish-start,' s'
      write(6,*)



!.....Residual normalization, convergence check  
      do i=1,nphi
        resor(i)=resor(i)*rnor(i)
      end do

      ! Write to monitor file - Old style
      !include 'simpleMonitorResiduals.h'


!.....Proveri kako stojimo sa rezidualom
      source=max(resor(iu),resor(iv),resor(iw),resor(ip),resor(ien),resor(icon)) 
      if(source.gt.slarge) then
          write(6,"(//,10x,a)") "*** Program terminated -  iterations diverge ***" 
          stop ! zavrsi program
      endif
!
!==========================================================
!.....False time steping: skoci na next line after the end of time loop
!==========================================================
      if(.not.ltransient) then
          if(source.lt.sormax) exit time_loop
      end if

      if(ltransient) then 

         ! Konverigao u okviru timestep-a ili potrosio sve iteracije:
         if(source.lt.sormax.or.iter.ge.maxit) then 

            if(const_mflux) then
              !# Correct driving force for a constant mass flow rate.
              include 'constant_mass_flow_forcing.f90'
            endif

            if(mod(itime,nzapis).eq.0) call writefiles
            call writehistory !<- write monitoring points
            call calc_statistics 


            ! # Create and save a frame for animation
            if(mod(itime,100).eq.0) then 
               include 'create_and_save_frame.f90'
            endif


            ! # Add fluctations by random perturbations to maintain turbulence in channel
            ! if(mod(itime,nzapis).eq.0) then
            !   call add_random_noise_to_field(u,5)
            !   call add_random_noise_to_field(v,5)
            !   call add_random_noise_to_field(w,5)
            !   write(6,'(a)') "Added random noise to velocity field!"
            ! endif

            cycle time_loop ! idi na pocetak 'time_loop'-a sa sledecom vrednosti indeksa
         
         endif

      end if 

      end do iteration_loop


!===============================================
!.....Write field values after nzapis iterations 
!===============================================
      if(.not.ltransient) then
        if(mod(itime,nzapis).eq.0) then
         call writefiles
         call write_restart_files
        endif
      endif

 
      if(ltransient) call flush(6)
 
      end do time_loop

!===============================================
!.....Write field values for the next run 
!===============================================
      call writefiles
      call write_restart_files
      

end program caffa3d

