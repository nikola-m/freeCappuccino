!
!***********************************************************************
!
  subroutine global_isum(i) 
!
!***********************************************************************
!
!   Estimates global sum among all processors of an integer variable.                    
!
!***********************************************************************
!
  use types
  
  implicit none

  include 'mpif.h'

  integer :: i

  integer :: isum
  integer :: ierr

  call mpi_allreduce      &               
   (i,                    & ! send buffer
    isum,                 & ! recv buffer 
    1,                    & ! length     
    mpi_integer,          & ! datatype  
    mpi_sum,              & ! operation 
    mpi_comm_world,       & ! communicator            
    ierr) 

  i = isum

  end subroutine global_isum