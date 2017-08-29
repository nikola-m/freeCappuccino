!
!***********************************************************************
!
  subroutine global_max(phi) 
!
!***********************************************************************
!
!   Estimates global maximum value of phi among all processes.
!   Every process creates it own maximum end stores it phi.
!   Than it all gets gathered and summed and sent back to each process.                       
!
!***********************************************************************
!
  use types
  
  implicit none

  include 'mpif.h'

  real(dp) :: phi

  real(dp) :: phimax
  integer :: ierr

  call mpi_allreduce      &               
   (phi,                  & ! send buffer
    phimax,               & ! recv buffer 
    1,                    & ! length     
    mpi_double_precision, & ! datatype  
    mpi_max,              & ! operation 
    mpi_comm_world,       & ! communicator            
    ierr) 

  phi = phimax

  end subroutine global_max