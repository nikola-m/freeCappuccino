!
!***********************************************************************
!
  subroutine abort_mission
!
!***********************************************************************
!
!   Finalize parallel communication and  program execution.                    
!
!***********************************************************************
!
  use parameters
  
  implicit none

  include 'mpif.h'

  call mpi_barrier(mpi_comm_world,ierr)
  call mpi_finalize(ierr)
  stop

  end subroutine