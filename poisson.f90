program poisson
!***********************************************************************
!
! Test discretisation of the Lapalcian operator on example problem:
! Laplacian(p) = -8 pi^2 sin(2 pi x) sin(2 pi y),
! with exact solution:
! p_exact = sin(2 pi x) sin(2 pi y)
! Test of meshes with different resolution and check the rate of grid 
! convergence.
!
!***********************************************************************
  use types
  use parameters
  use geometry
  use sparse_matrix, only: create_CSR_matrix_from_mesh_data,su,sv
  use utils, only: show_logo
  use LIS_linear_solver_library

  implicit none
!
!***********************************************************************
!

! 
! Local variables 
!
  integer :: i, ijb, iface, ijp, ijn
  real(dp) :: lh
  real(dp), dimension(:), allocatable:: p

!
!***********************************************************************
!

!
! > Open log file
! 
  open(unit=6,file='log_poisson')
  rewind 6
! 
! > Print code logo and timestamp in log file
!
  call show_logo

!
! >  Open & Read mesh file, calculate mesh geometrical quantities, allocate arrays
!
  call mesh_geometry

  allocate(p(numTotal))

!
! >  Index arrays of matrix elements stored in CSR format
!
  call create_CSR_matrix_from_mesh_data

!
! > Discretisation and solution of the tet problem
!

  ! Source term
  su = 0.0_dp
  su(1:numCells) = 8*pi**2*sin(2*pi*xc(1:numCells))*sin(2*pi*yc(1:numCells))*Vol(1:numCells)

  ! Initialize solution
  p(1:numTotal) = 0.0_dp

  do i=1,numBoundaryFaces
  iface = numInnerFaces+i
  ijb = numCells+i
  p(ijb) = 0.0_dp 
  enddo

  !  Coefficient array for Laplacian
  sv = -1.0_dp       

  ! Laplacian operator and BCs         
  call fvm_laplacian(sv,p) 

  sor(ip) = 1e-16
  nsw(ip) = 1000

  ! Solve system
  write(*,'(a)') ' '
  ! 1)
  ! call iccg(p,ip) 
  ! 2)
  ! call gaussSeidel(p,ip)
  ! 3)
  call solve_csr(numCells,nnz,ioffset,ja,a,su,p) 
  write(*,'(a)') ' '

  do i=1,numCells
    write(6,'(es11.4)') p(i)
  enddo 

  ! Cell size
  ijp = owner(1)
  ijn  = neighbour(1)
  lh = abs(xc(ijp)-xc(ijn))
  
  ! L_inf error norm: 
  write(*,'(a)') ' '
  write(*,'(2es11.4)') lh , maxval( abs( p(1:numCells)-sin(2*pi*xc(1:numCells))*sin(2*pi*yc(1:numCells)) ) )


end program