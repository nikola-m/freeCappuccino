!***********************************************************************
!
subroutine linear_solver(fi,ifi)
!
!***********************************************************************
!
!
!
  use types
  use parameters
  use geometry
  use sparse_matrix
  use title_mod

  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: ifi
  real(dp), dimension(numTotal), intent(inout) :: fi 
!
! Local variables
!

  if( linSolver(ifi) .eq. xx ) then

    call dpcg(fi, ifi)

  elseif( linSolver(ifi) .eq. xx ) then 

    call iccg(fi, ifi)

  elseif( linSolver(ifi) .eq. xx ) then  

    call bicgstab(fi, ifi) 

  elseif( linSolver(ifi) .eq. xx ) then 

    call solve_csr( numCells, nnz, ioffset, ja , a, su, fi)

  elseif( linSolver(ifi) .eq. xx ) then 

    call pmgmres_ilu ( numCells, nnz, ioffset, ja, a, diag, fi(1:numCells), ip, su, 100, 4, 1e-8, sor(ip) )
    
  endif


end subroutine