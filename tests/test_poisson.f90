program poisson
!***********************************************************************
!
! Test discretisation of the Lapalcian operator on example problem:
! -Lapalcian(p) = -8 pi^2 sin(2 pi x) sin(2 pi y),
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

  implicit none
!
!***********************************************************************
!

! 
! Local variables 
!
  integer :: i, ijp, ijn, inp, iface
  integer :: nsw_backup
  real(dp) :: sor_backup

!
!***********************************************************************
!

!
! > Open log file
! 
  open(unit=6,file='log_test_poisson.txt')
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
  su(1:numCells) = -8*pi**2*sin(2*pi*xc(1:numCells))*sin(2*pi*yc(1:numCells))*Vol(1:numCells)

  ! Initialize solution
  p(1:numCells) = 0.0_dp

  do i=1,numBoundaryFaces
  iface = numInnerFaces+i
  ijb = numCells+i
  p(ijb) = sin(2*pi*xf(iface))*sin(2*pi*yf(iface))
  enddo

  !  Coefficient array for Laplacian
  sv = 1.0_dp       

  ! Laplacian operator and BCs         
  call fvm_laplacian(sv,p) 

  sor_backup = sor(ip)
  nsw_backup = nsw(ip)

  sor(ip) = 1e-16
  nsw(ip) = 1000

  ! Solve system
  write(*,'(a)') ' '
  call iccg(p,ip) 
  ! call gaussSeidel(p,ip) 

  do i=1,numCells
    write(6,'(es11.4)') p(i)
  enddo 

  write(*,'(a)') ' '
  write(*,'(a,es11.4)') 'L_inf error norm: ', maxval( abs( p(1:numCells)-sin(2*pi*xc(1:numCells))*sin(2*pi*yc(1:numCells)) ) )


end program