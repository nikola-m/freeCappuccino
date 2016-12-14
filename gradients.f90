module gradients

use types
use parameters
use geometry, only:numTotal,numCells

implicit none

logical :: lstsq, lstsq_qr, lstsq_dm, gauss     ! gradient discretization approach

real(dp),dimension(:,:), allocatable ::  dmat  !  d(6,nxyz) - when using bn, or dm version of the subroutine
real(dp),dimension(:,:,:), allocatable ::  dmatqr  ! when using qr version of the subroutine size(3,6,nxyz)!


interface grad
  module procedure grad_scalar_field
  module procedure grad_vector_field
end interface


private

public :: lstsq, lstsq_qr, lstsq_dm, gauss
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
  real(dp), dimension(3,numTotal), intent(inout) :: dPhidxi

  dPhidxi(:,:) = 0.0_dp

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
  real(dp), dimension(3,numTotal), intent(inout) :: dPhidxi

  dPhidxi(:,:) = 0.0_dp
  
  if (lstsq) then 
    call grad_lsq(phi,dPhidxi,2,dmat)
  elseif (lstsq_qr) then 
    call grad_lsq_qr(phi,dPhidxi,2,dmatqr)
  elseif (lstsq_dm) then 
    call grad_lsq_dm(phi,dPhidxi,2,dmat)
  else
    call grad_gauss(phi,dPhidxi(1,:),dPhidxi(2,:),dPhidxi(3,:))
  endif 

end subroutine


subroutine grad_vector_field(U,V,W,dUdxi,dVdxi,dWdxi)

use types
use parameters

implicit none

  real(dp), dimension(numTotal), intent(in) :: U,V,W
  real(dp), dimension(3,numTotal), intent(inout) :: dUdxi,dVdxi,dWdxi

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



! Least square gradients
include 'grad_lsq.f90'


! Least square gradients via QR decomposition
include 'grad_lsq_qr.f90'


! Weighted least square gradients
include 'grad_lsq_dm.f90'


! Gauss gradients
include 'grad_gauss.f90'


end module gradients