module gradients

use types
use parameters
use geometry, only:numTotal,numCells

implicit none

logical :: lstsq, lstsq_qr, lstsq_dm, gauss       ! Gradient discretization approach
character(len=20) :: limiter = 'none' ! Gradient limiter. Options: none, Barth-Jespersen, Venkatakrishnan, MVenkatakrishnan

real(dp),dimension(:,:), allocatable ::  dmat  !  d(6,nxyz) - when using bn, or dm version of the subroutine
real(dp),dimension(:,:,:), allocatable ::  dmatqr  ! when using qr version of the subroutine size(3,6,nxyz)!


interface grad
  module procedure grad_scalar_field
  module procedure grad_vector_field
end interface


private

public :: lstsq, lstsq_qr, lstsq_dm, gauss
public :: limiter
public :: grad,allocate_gradients, &
          create_lsq_gradients_matrix



contains


subroutine allocate_gradients
  implicit none
  
  integer :: ierr

  if( lstsq .or. lstsq_dm ) then
    allocate(dmat(6,numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: dmat"
  elseif( lstsq_qr ) then
    allocate(dmatqr(3,6,numCells),stat=ierr) 
      if(ierr /= 0)write(*,*)"allocation error: dmatqr"
  endif


end subroutine



subroutine create_lsq_gradients_matrix(phi,dPhidxi)
!
!  Discussion:
!    Prepare System Matrix For Least-Squares Gradient Calculation.
!    It is done by setting this --v value to one.
!           call grad_lsq(U,dUdxi,1,D)
!
use types
use parameters

implicit none

  real(dp), dimension(numTotal), intent(in) :: phi
  real(dp), dimension(3,numCells), intent(inout) :: dPhidxi

  if (lstsq) then
    call grad_lsq(phi,dPhidxi,1,dmat)
  elseif (lstsq_qr) then 
    call grad_lsq_qr(phi,dPhidxi,1,dmatqr)
  elseif (lstsq_dm) then 
    call grad_lsq_dm(phi,dPhidxi,1,dmat)
  endif 

end subroutine



subroutine grad_scalar_field(phi,dPhidxi)

use types
use parameters

implicit none

  real(dp), dimension(numTotal), intent(in) :: phi
  real(dp), dimension(3,numCells), intent(inout) :: dPhidxi

  dPhidxi = 0.0_dp
  
  if (lstsq) then 
    call grad_lsq(phi,dPhidxi,2,dmat)
  elseif (lstsq_qr) then 
    call grad_lsq_qr(phi,dPhidxi,2,dmatqr)
  elseif (lstsq_dm) then 
    call grad_lsq_dm(phi,dPhidxi,2,dmat)
  else
    call grad_gauss(phi,dPhidxi(1,:),dPhidxi(2,:),dPhidxi(3,:))
  endif 

  ! Gradient limiter:
  if(adjustl(limiter) == 'Barth-Jespersen') then
    call slope_limiter_Barth_Jespersen(phi, dPhidxi)
  elseif(adjustl(limiter) == 'Venkatakrishnan') then
    call slope_limiter_Venkatakrishnan(phi, dPhidxi)
  elseif(adjustl(limiter) == 'MVenkatakrishnan') then
    call slope_limiter_modified_Venkatakrishnan(phi, dPhidxi)
  else
    ! no-limit
  endif

end subroutine


subroutine grad_vector_field(U,V,W,dUdxi,dVdxi,dWdxi)

use types
use parameters

implicit none

  real(dp), dimension(numTotal), intent(in) :: U,V,W
  real(dp), dimension(3,numCells), intent(inout) :: dUdxi,dVdxi,dWdxi

  dUdxi=0.0_dp
  dVdxi=0.0_dp
  dWdxi=0.0_dp

  if (lstsq) then
    call grad_lsq(U,dUdxi,2,dmat)
    call grad_lsq(V,dVdxi,2,dmat)
    call grad_lsq(W,dWdxi,2,dmat)
  elseif (lstsq_qr) then
    call grad_lsq_qr(U,dUdxi,2,dmatqr)
    call grad_lsq_qr(V,dVdxi,2,dmatqr)
    call grad_lsq_qr(W,dWdxi,2,dmatqr)
  elseif (lstsq_dm) then
    call grad_lsq_dm(U,dUdxi,2,dmat)
    call grad_lsq_dm(V,dVdxi,2,dmat)
    call grad_lsq_dm(W,dWdxi,2,dmat)
  else
    call grad_gauss(U,dUdxi(1,:),dUdxi(2,:),dUdxi(3,:))
    call grad_gauss(V,dVdxi(1,:),dVdxi(2,:),dVdxi(3,:))
    call grad_gauss(W,dWdxi(1,:),dWdxi(2,:),dWdxi(3,:))
  endif

end subroutine


!***********************************************************************
!
subroutine slope_limiter_modified_Venkatakrishnan(phi, dPhidxi)
!
!***********************************************************************
!
!     Calculates slope limiter and appiies to scalar gradient:
!     Wang modified Venkatakrishnan slope limiter
!     Ref.: Z. J. Wang. "A Fast Nested Multi-grid Viscous Flow Solver for Adaptive Cartesian/Quad Grids",
!     International Journal for Numerical Methods in Fluids. 33. 657–680. 2000.
!     The same slope limiter is used in Fluent.
!
!***********************************************************************
!
  use types
  use parameters
  use sparse_matrix, only: ioffset,ja,diag
  use geometry, only: numTotal,numCells,xc,yc,zc

  implicit none

  !     Input
  real(dp),dimension(numTotal) :: phi
  real(dp),dimension(3,numCells) :: dPhidxi


  !     Locals
  integer :: inp,ijp,ijn,k
  real(dp) :: phi_p
  real(dp) :: cell_neighbour_value,gradfiXdr,slopelimit
  real(dp) :: deltam,deltap,epsi
  real(dp) :: phi_max,phi_min
  real(dp) :: fimax,fimin


  fimin = minval(phi(1:numCells))
  fimax = maxval(phi(1:numCells))

  do inp = 1, numCells

    ! Values at cell center:
    phi_p = phi(inp)

    ! max and min values over current cell and neighbors
    phi_max = phi(ja( ioffset(inp) ))
    phi_min = phi(ja( ioffset(inp) ))

    do k=ioffset(inp)+1, ioffset(inp+1)-1
      phi_max = max( phi_max, phi(ja(k)) )
      phi_min = min( phi_max, phi(ja(k)) )      
    enddo


    slopelimit = 1.0_dp

    do k=ioffset(inp), ioffset(inp+1)-1

      if (k == diag(inp)) cycle
   
      ijp = inp
      ijn = ja(k)

      gradfiXdr=dPhidxi(1,ijp)*(xc(ijn)-xc(ijp))+dPhidxi(2,ijp)*(yc(ijn)-yc(ijp))+dPhidxi(3,ijp)*(zc(ijn)-zc(ijp)) 

      ! Find unlimited value:
      cell_neighbour_value =  phi_p + gradfiXdr 


      deltam = cell_neighbour_value - phi_p
      if (deltam .gt. 0.0d0) then
          deltap = phi_max-phi_p
      else
          deltap = phi_min-phi_p
      endif

      ! Wang proposition for epsilon
      epsi = (0.05*( fimax-fimin ))**2 

      slopelimit = max(min(slopelimit, 1./(deltam+small)*((deltap+epsi)*deltam+2*deltam**2*deltap) &
                                                     /(deltap**2+2*deltam**2+deltap*deltam+epsi+small)),zero)

    enddo
    !print*,slopelimit
    dPhidxi(:,inp) = slopelimit*dPhidxi(:,inp)

  enddo

end subroutine



!***********************************************************************
!
subroutine slope_limiter_Barth_Jespersen(phi, dPhidxi)
!
!***********************************************************************
!
!     Calculates slope limiter and appiies to scalar gradient:
!     Barth and Jespersen slope limiter:
!
!     AIAA-89-0366, The design and application of upwind schemes
!     on unstructured meshes, T.J.Barth, D.C.Jespersen, 1989.
!
!***********************************************************************
!
  use types
  use parameters
  use sparse_matrix, only: ioffset,ja,diag
  use geometry, only: numTotal,numCells,xc,yc,zc

  implicit none

  !     Input
  real(dp),dimension(numTotal) :: phi
  real(dp),dimension(3,numCells) :: dPhidxi


  !     Locals
  integer :: inp,ijp,ijn,k
  real(dp) :: phi_p
  real(dp) :: slopelimit
  real(dp) :: delta_face
  real(dp) :: phi_max,phi_min,r
  real(dp) :: fimax,fimin,deltamax,deltamin


  fimin = minval(phi(1:numCells))
  fimax = maxval(phi(1:numCells))

  do inp = 1, numCells

    ! Values at cell center:
    phi_p = phi(inp)

    ! max and min values over current cell and neighbors
    phi_max = phi(ja( ioffset(inp) ))
    phi_min = phi(ja( ioffset(inp) ))

    do k=ioffset(inp)+1, ioffset(inp+1)-1
      phi_max = max( phi_max, phi(ja(k)) )
      phi_min = min( phi_max, phi(ja(k)) )      
    enddo


    deltamax = fimax - phi(inp)
    deltamin = fimin - phi(inp)

    slopelimit = 1.0_dp

    do k=ioffset(inp), ioffset(inp+1)-1

      if (k == diag(inp)) cycle
   
      ijp = inp
      ijn = ja(k)

      delta_face=dPhidxi(1,ijp)*(xc(ijn)-xc(ijp))+dPhidxi(2,ijp)*(yc(ijn)-yc(ijp))+dPhidxi(3,ijp)*(zc(ijn)-zc(ijp)) 


     if( abs(delta_face) < 1.e-6 )then
       r = 1.0_dp
     else if( delta_face > 0.0 )then
       r = deltamax/delta_face
     else
       r = deltamin/delta_face
     endif

     slopelimit = min( slopelimit , r )

    enddo
    !print*,slopelimit
    dPhidxi(:,inp) = slopelimit*dPhidxi(:,inp)

  enddo

end subroutine




!***********************************************************************
!
subroutine slope_limiter_Venkatakrishnan(phi, dPhidxi)
!
!***********************************************************************
!
!     Calculates slope limiter and appiies to scalar gradient:
!     Venkatakrishnan slope limiter:
!
!    AIAA-93-0880, On the accuracy of limiters and convergence
!    to steady state solutions, V.Venkatakrishnan, 1993
!
!***********************************************************************
!
  use types
  use parameters
  use sparse_matrix, only: ioffset,ja,diag
  use geometry, only: numTotal,numCells,xc,yc,zc

  implicit none

  !     Input
  real(dp),dimension(numTotal) :: phi
  real(dp),dimension(3,numCells) :: dPhidxi


  !     Locals
  integer :: inp,ijp,ijn,k
  real(dp) :: phi_p
  real(dp) :: slopelimit
  real(dp) :: delta_face
  real(dp) :: phi_max,phi_min,r
  real(dp) :: fimax,fimin,deltamax,deltamin


  fimin = minval(phi(1:numCells))
  fimax = maxval(phi(1:numCells))

  do inp = 1, numCells

    ! Values at cell center:
    phi_p = phi(inp)

    ! max and min values over current cell and neighbors
    phi_max = phi(ja( ioffset(inp) ))
    phi_min = phi(ja( ioffset(inp) ))

    do k=ioffset(inp)+1, ioffset(inp+1)-1
      phi_max = max( phi_max, phi(ja(k)) )
      phi_min = min( phi_max, phi(ja(k)) )      
    enddo


    deltamax = fimax - phi(inp)
    deltamin = fimin - phi(inp)

    slopelimit = 1.0_dp

    do k=ioffset(inp), ioffset(inp+1)-1

      if (k == diag(inp)) cycle
   
      ijp = inp
      ijn = ja(k)

      delta_face=dPhidxi(1,ijp)*(xc(ijn)-xc(ijp))+dPhidxi(2,ijp)*(yc(ijn)-yc(ijp))+dPhidxi(3,ijp)*(zc(ijn)-zc(ijp)) 


     if( abs(delta_face) < 1.e-6 )then
       r = 1.0_dp
     else if( delta_face > 0.0 )then
       r = deltamax/delta_face
     else
       r = deltamin/delta_face
     endif

     slopelimit = min( slopelimit , (r**2+2.0*r)/(r**2+r+2.0) )

    enddo

    dPhidxi(:,inp) = slopelimit*dPhidxi(:,inp)

  enddo

end subroutine



! Least square gradients
include 'grad_lsq.f90'


! Least square gradients via QR decomposition
include 'grad_lsq_qr.f90'


! Weighted least square gradients
include 'grad_lsq_dm.f90'


! Gauss gradients
include 'grad_gauss.f90'


end module gradients