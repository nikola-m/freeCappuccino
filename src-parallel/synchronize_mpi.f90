!
!***********************************************************************
!
  subroutine synchronize 
!
!***********************************************************************
!
!  Calls mpi_barrier and syncronizes all processes at an end of a task.                      
!
!***********************************************************************
!
  use types
  
  implicit none

  include 'mpif.h'

  integer :: ierr

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

  end subroutine