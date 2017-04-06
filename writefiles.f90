!***********************************************************************
!
subroutine writefiles
!
!***********************************************************************
!
  use types
  use parameters
  use title_mod
  use geometry
  use variables
  use statistics
  use sparse_matrix

  implicit none
!
!***********************************************************************
!
  integer :: i
  integer :: output_unit
  character(len=6) :: timechar

  ! Write in a char variable current timestep number and create a folder with this name
  write(timechar,'(i6)') itime
  call execute_command_line("mkdir "//trim(out_folder_path)//"/"//adjustl(timechar))



  ! Open an output file in this folder for velocity
  call get_unit ( output_unit )
  open ( unit = output_unit, file = trim(out_folder_path)//'/'//trim(adjustl(timechar))//'/U')

  rewind output_unit

  ! Write output to file
  do i=1,numCells
    write(output_unit,'(a,es11.4,1x,es11.4,1x,es11.4,a)') '(',u(i),v(i),w(i),')'
  enddo

  close(output_unit)



  ! Open an output file in this folder for pressure
  call get_unit ( output_unit )
  open ( unit = output_unit, file = trim(out_folder_path)//'/'//trim(adjustl(timechar))//'/p')

  rewind output_unit

  ! Write output to file
  do i=1,numCells
    write(output_unit ,'(es11.4)') p(i)
  enddo

  close(output_unit)



  if(solveTKE) then

    ! Open an output file in this folder for Turbulence kinetic energy
    call get_unit ( output_unit )
    open ( unit = output_unit, file = trim(out_folder_path)//'/'//trim(adjustl(timechar))//'/k')

    rewind output_unit

    ! Write output to file
    do i=1,numCells
      write(output_unit ,'(es11.4)') te(i)
    enddo

    close(output_unit)

  endif



  if( solveEpsilon ) then

    ! Open an output file in this folder for Epsilon
    call get_unit ( output_unit )
    open ( unit = output_unit, file = trim(out_folder_path)//'/'//trim(adjustl(timechar))//'/epsilon')

    rewind output_unit

    ! Write output to file
    do i=1,numCells
      write(output_unit ,'(es11.4)') ed(i)
    enddo

    close(output_unit)

  endif

  

  if( solveOmega ) then

    ! Open an output file in this folder for Omega
    call get_unit ( output_unit )
    open ( unit = output_unit, file = trim(out_folder_path)//'/'//trim(adjustl(timechar))//'/omega')

    rewind output_unit

    ! Write output to file
    do i=1,numCells
      write(output_unit ,'(es11.4)') ed(i)
    enddo

    close(output_unit)

  endif      


! !----------------------------------------------------------
! !     [TECPLOT FILE FORMAT: ]
! !----------------------------------------------------------
!       OPEN(UNIT=82,FILE=trim(out_folder_path)//'/TECPLOT_FILE.plt')
!       REWIND 82

!       WRITE(82,*) 'TITLE     = " "'
!       WRITE(82,*) 'VARIABLES = "X"'
!       WRITE(82,*) '"Y"'
!       WRITE(82,*) '"Z"'
!       WRITE(82,*) '"U"'
!       WRITE(82,*) '"V"'
!       WRITE(82,*) '"W"'
!       WRITE(82,*) '"P"'
!       WRITE(82,*) '"TE"'
!       WRITE(82,*) '"ED"'
!       WRITE(82,*) '"VIS"'
!       WRITE(82,*) '"UU"'
!       WRITE(82,*) '"VV"'
!       WRITE(82,*) '"WW"'
!       WRITE(82,*) '"Q"'
!       WRITE(82,*) 'ZONE T=" "'
!       WRITE(82,*) 'I=',NI, ' ,J=',NJ, ' ,K=',NK,', F=POINT'

!       DO K=1,NK
!       DO J=1,NJ
!       DO I=1,NI
!       IJK=(K-1)*NI*NJ+(I-1)*NJ+J
!       INBC=(I-1)*NJ+J 
!       QVortex = 0.5*(Vorticity(ijk)-magStrain(ijk)) ! Q criteria for vortex identification
!       WRITE(82,*) XC(IJK),YC(IJK),ZC(IJK),U(IJK),V(IJK),W(IJK), &
!                  P(IJK),TE(IJK),ED(IJK),VIS(IJK), &
!                  UU(IJK),VV(IJK),WW(IJK),QVortex
!       END DO
!       END DO
!       END DO

!       REWIND 82
!       CLOSE (82)

!       if (ltransient) then
! !--------------------------------------------------------------
! !    [ writing of the statistics tecplot file]
! !----------------------------------------------------------
!       OPEN(UNIT=87,FILE=trim(out_folder_path)//'/TECPLOT_STAT.plt')                 ! <- statistics postprocessing file
!       REWIND 87

!       WRITE(87,*) 'TITLE     = " "'
!       WRITE(87,*) 'VARIABLES = "X"'
!       WRITE(87,*) '"Y"'
!       WRITE(87,*) '"Z"'
!       WRITE(87,*) '"U_AVER"'
!       WRITE(87,*) '"V_AVER"'
!       WRITE(87,*) '"W_AVER"'
!       WRITE(87,*) '"UU_AVER"'
!       WRITE(87,*) '"VV_AVER"'
!       WRITE(87,*) '"WW_AVER"'
!       WRITE(87,*) '"UV_AVER"'
!       WRITE(87,*) '"UW_AVER"'
!       WRITE(87,*) '"VW_AVER"'
!       WRITE(87,*) '"TE_AVER"'
!       WRITE(87,*) 'ZONE T=" "'
!       WRITE(87,*) 'I=',NI, ' ,J=',NJ, ' ,K=',NK,', F=POINT'
!       DO K=1,NK
!       DO J=1,NJ
!       DO I=1,NI
!       IJK=(K-1)*NI*NJ+(I-1)*NJ+J
!       WRITE(87,*) XC(IJK),YC(IJK),ZC(IJK), &
!                   U_AVER(IJK),V_AVER(IJK),W_AVER(IJK), &   
!                   UU_AVER(IJK),VV_AVER(IJK),WW_AVER(IJK), &
!                   UV_AVER(IJK),UW_AVER(IJK),VW_AVER(IJK), &
!                   TE_AVER(IJK)
!       END DO
!       END DO
!       END DO

!       CLOSE(87)
! !--------------------------------------------------------------
!       endif

end subroutine
