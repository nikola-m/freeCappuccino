!
!***********************************************************************
!
  subroutine global_min(phi) 
!
!***********************************************************************
!
!   Estimates global minimum value of phi among all processes.
!   Every process creates it own minimum end stores it phi.
!   Than it all gets gathered and summed and sent back to each process.                       
!
!***********************************************************************
!
  use types
  
  implicit none

  include 'mpif.h'

  real(dp) :: phi

  real(dp) :: phimin
  integer :: ierr

  call mpi_allreduce      &               
   (phi,                  & ! send buffer
    phimin,               & ! recv buffer 
    1,                    & ! length     
    mpi_double_precision, & ! datatype  
    mpi_min,              & ! operation 
    mpi_comm_world,       & ! communicator            
    ierr) 

  phi = phimin

  end subroutine global_min