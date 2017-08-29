!
!***********************************************************************
!
  subroutine global_sum(phi) 
!
!***********************************************************************
!
!   Estimates global sum among all processors.
!   Used e.g. in dot product. Every process creates it own sum
!   end stores it in phi. Than it all gets gathered and summed and sent
!   back to each process.                       
!
!***********************************************************************
!
  use types
  
  implicit none

  include 'mpif.h'

  real(dp) :: phi

  real(dp) :: phisum
  integer :: ierr

  call mpi_allreduce      &               
   (phi,                  & ! send buffer
    phisum,               & ! recv buffer 
    1,                    & ! length     
    mpi_double_precision, & ! datatype  
    mpi_sum,              & ! operation 
    mpi_comm_world,       & ! communicator            
    ierr) 

  phi = phisum

  end subroutine global_sum