module gradients
!
! Module for cell center gradients, gradient limiters and surface normal gradients.
!

use types
use parameters
use geometry, only: numCells,numPCells,numTotal,npro,xc,yc,zc
use sparse_matrix, only: ioffset,ja,diag

implicit none

logical :: lstsq, lstsq_qr, lstsq_dm, gauss        ! Gradient discretization approach
character(len=20) :: limiter                       ! Gradient limiter. Options: none, Barth-Jespersen, Venkatakrishnan, mVenkatakrishnan
character(len=12) :: sngrad_corr                   ! Surface normal gradient correction scheme. Options: Skewness, Offset, Uncorrected.

real(dp),dimension(:,:), allocatable ::  dmat      !  d(6,nxyz) - when using bn, or dm version of the subroutine
real(dp),dimension(:,:,:), allocatable ::  dmatqr  !  when using qr version of the subroutine size(3,6,nxyz)!


interface grad
  module procedure grad_scalar_field
  module procedure grad_vector_field
end interface

interface sngrad
  module procedure sngrad_scalar_field
  module procedure sngrad_vector_field
end interface


private

public :: lstsq, lstsq_qr, lstsq_dm, gauss,limiter,sngrad_corr
public :: grad,sngrad,allocate_gradients,create_lsq_gradients_matrix



contains


!***********************************************************************
!
subroutine allocate_gradients
!
!***********************************************************************
!
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


!***********************************************************************
!
subroutine create_lsq_gradients_matrix(phi,dPhidxi)
!
!***********************************************************************
!
!  Discussion:
!    Prepare System Matrix For Least-Squares Gradient Calculation.
!    It is done by setting this --v value to one.
!           call grad_lsq(U,dUdxi,1,D)
!
!***********************************************************************
!

implicit none

  real(dp), dimension(numTotal), intent(in) :: phi
  real(dp), dimension(3,numPCells), intent(inout) :: dPhidxi

  if (lstsq) then
    call grad_lsq(phi,dPhidxi,1,dmat)
  elseif (lstsq_qr) then 
    call grad_lsq_qr(phi,dPhidxi,1,dmatqr)
  elseif (lstsq_dm) then 
    call grad_lsq_dm(phi,dPhidxi,1,dmat)
  endif 

end subroutine


!***********************************************************************
!
subroutine grad_scalar_field(phi,dPhidxi)
!
!***********************************************************************
!

implicit none

  real(dp), dimension(numTotal), intent(inout) :: phi
  real(dp), dimension(3,numPCells), intent(inout) :: dPhidxi

  ! Before the calculation of the gradients we exchange the information
  ! between processes to make sure that the freshest field values are in the buffer.
  ! MPI exchange:
  call exchange(phi)

  ! Initialize gradient:
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
  if( trim(adjustl(limiter)) == 'Barth-Jespersen') then

    call slope_limiter_Barth_Jespersen( phi, dPhidxi )

  elseif( trim(adjustl(limiter)) == 'Venkatakrishnan') then

    call slope_limiter_Venkatakrishnan( phi, dPhidxi )

  elseif( trim(adjustl(limiter)) == 'mVenkatakrishnan') then

    call slope_limiter_modified_Venkatakrishnan( phi, dPhidxi )
    
  else
    ! no-limit
  endif

  ! MPI exchange:
  call exchange( dPhidxi(1,:) )
  call exchange( dPhidxi(2,:) )
  call exchange( dPhidxi(3,:) )

end subroutine


!***********************************************************************
!
subroutine grad_vector_field(U,V,W,dUdxi,dVdxi,dWdxi)
!
!***********************************************************************
!

  implicit none

  real(dp), dimension(numTotal), intent(inout) :: U,V,W
  real(dp), dimension(3,numPCells), intent(inout) :: dUdxi,dVdxi,dWdxi

  ! Before the calculation of the gradients we exchange the information
  ! between processes to make sure that the freshest field values are in the buffer.
  ! MPI exchange:
  call exchange(U)
  call exchange(V)
  call exchange(W)

  ! Initialize gradient:
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


  ! MPI exchange:
  call exchange_short( dUdxi(1,:) )
  call exchange_short( dUdxi(2,:) )
  call exchange_short( dUdxi(3,:) )

  call exchange_short( dVdxi(1,:) )
  call exchange_short( dVdxi(2,:) )
  call exchange_short( dVdxi(3,:) )

  call exchange_short( dWdxi(1,:) )
  call exchange_short( dWdxi(2,:) )
  call exchange_short( dWdxi(3,:) )

end subroutine


!***********************************************************************
!
subroutine set_phi_min_max(phi)
!
!***********************************************************************
!
! Calculates and stores for every cell, a minimum and maximum value
! of field variable PHI, over current cell and its neighbours.
! PHI_MAX and PHI_MIN are used for gradient limiter calculation.
!
! NOTE: Buffer needs to have fresh values! 
! Exchange phi should be done previously.
!
!***********************************************************************
!

  use geometry, only: numCells,numInnerFaces,numTotal,npro,iProcFacesStart,iProcStart,owner,neighbour
  use variables, only: phimax,phimin

  implicit none

  ! Input
  real(dp),dimension(numTotal) :: phi

  ! Locals
  integer :: i,ijp,ijn,iface

  ! Initialize max and min arrays with values of phi in each cell
  phimax = phi(1:numCells)
  phimin = phi(1:numCells)

  ! > Loop over neighbours 

  ! Inner faces
  do i = 1,numInnerFaces

    ijp = owner(i)
    ijn = neighbour(i)

    phimax(ijp) = max( phimax(ijp), phi(ijn) )
    phimax(ijn) = max( phimax(ijn), phi(ijp) )


    phimin(ijp) = min( phimin(ijp), phi(ijn) )
    phimin(ijn) = min( phimin(ijn), phi(ijp) )

  enddo

  ! Faces on processor boundary
  do i=1,npro

    iface = iProcFacesStart + i
    ijp = owner( iface )
    ijn = iProcStart + i

    phimax(ijp) = max( phimax(ijp), phi(ijn) )

    phimin(ijp) = min( phimin(ijp), phi(ijn) )

  enddo

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
  use variables, only:phimax,phimin

  implicit none

  !     Input
  real(dp),dimension(numTotal) :: phi
  real(dp),dimension(3,numPCells) :: dPhidxi


  !     Locals
  integer :: inp,ijp,ijn,k

  ! Look at the reference epsprim \in [0.01,0.2]
  real(dp), parameter :: epsprim = 0.05_dp

  real(dp) :: phi_p
  real(dp) :: cell_neighbour_value,gradfiXdr,slopelimit
  real(dp) :: deltam,deltap,epsi
  real(dp) :: glomax,glomin


  ! Set global minimum and maximum; requires mpi_reduce
  glomin = minval(phi(1:numCells))
  glomax = maxval(phi(1:numCells))

  call global_min(glomin)
  call global_max(glomax)

  ! Sets phimin and phi max which are min and max value in the range of a cell and its neighbours.  
  call set_phi_min_max( phi )


  do inp = 1, numCells

    ! Values at cell center:
    phi_p = phi(inp)


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
          deltap = phimax(inp)-phi_p
      else
          deltap = phimin(inp)-phi_p
      endif

      ! Wang proposition for epsilon
      epsi =  epsprim*( glomax-glomin ) 
      slopelimit = max(                                                                          &
                        min(                                                                     &
                              slopelimit,                                                        &
                              1./(deltam+small)*((deltap**2+epsi**2)*deltam+2*deltam**2*deltap)  &
                                            /(deltap**2+2*deltam**2+deltap*deltam+epsi**2+small) &
                            ),                                                                   &
                        zero                                                                     &
                      )


    enddo

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

  implicit none

  !     Input
  real(dp),dimension(numTotal) :: phi
  real(dp),dimension(3,numPCells) :: dPhidxi


  !     Locals
  integer :: inp,ijp,ijn,k
  integer :: istart,iend
  real(dp) :: phi_p
  real(dp) :: slopelimit
  real(dp) :: delta_face
  real(dp) :: phi_max,phi_min,r
  real(dp) :: fimax,fimin,deltamax,deltamin


  fimin = minval(phi(1:numCells))
  fimax = maxval(phi(1:numCells))

  call global_min(fimin)
  call global_max(fimax)


  do inp = 1, numCells

    ! Values at cell center:
    phi_p = phi(inp)

    ! max and min values over current cell and neighbors
    istart = ioffset(inp)
    iend = ioffset(inp+1)-1

    phi_max = phi( ja( istart ) )
    phi_min = phi( ja( istart ) )

    do k = istart+1, iend      
      phi_max = max( phi_max, phi( ja(k) ) )
      phi_min = min( phi_max, phi( ja(k) ) )      
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

  implicit none

  !     Input
  real(dp),dimension(numTotal) :: phi
  real(dp),dimension(3,numPCells) :: dPhidxi


  !     Locals
  integer :: inp,ijp,ijn,k
  integer :: istart, iend
  real(dp) :: phi_p
  real(dp) :: slopelimit
  real(dp) :: delta_face
  real(dp) :: phi_max,phi_min,r
  real(dp) :: fimax,fimin,deltamax,deltamin


  fimin = minval(phi(1:numCells))
  fimax = maxval(phi(1:numCells))

  call global_min(fimin)
  call global_max(fimax)

  
  do inp = 1, numCells

    ! Values at cell center:
    phi_p = phi(inp)

    ! max and min values over current cell and neighbors
    istart = ioffset(inp)
    iend = ioffset(inp+1)-1

    phi_max = phi( ja( istart ) )
    phi_min = phi( ja( istart ) )

    do k = istart+1, iend      
      phi_max = max( phi_max, phi( ja(k) ) )
      phi_min = min( phi_max, phi( ja(k) ) )      
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


!******************************************************************************
!
subroutine sngrad_scalar_field(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
                               Fi, dFidxi, nrelax, approach, dfixi, dfiyi, dfizi, &
                               dfixii, dfiyii, dfizii)
!
!******************************************************************************
!
!  Surface normal gradient with non-orthogonal correction done in two
!  possible ways - either by skewness correction of intersection point
!  offset.
!
!  Check out reference paper: 
!    Mirkov, Rasuo, Kenjeres, JCP, Vol. 287, 2015.
!
!******************************************************************************
!
  implicit none
!
!******************************************************************************
! 
  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numPCells), intent(in) :: dFidxi
  integer, intent(in) :: nrelax
  character(len=12) :: approach
  real(dp), intent(out) :: dfixi, dfiyi, dfizi, dfixii, dfiyii, dfizii
!
! Locals
!
  real(dp) :: are,vole
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz
  real(dp) :: ixi1,ixi2,ixi3
  real(dp) :: dpn,costheta,costn

  real(dp) :: d1x,d1y,d1z
  real(dp) :: d2x,d2y,d2z
  real(dp) :: fxp,fxn

  real(dp) :: xpp,ypp,zpp,xep,yep,zep,xpnp,ypnp,zpnp,volep
  real(dp) :: nablaFIxdnnp,nablaFIxdppp

!
!******************************************************************************
!

  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0_dp-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! Components of the unit vector i_ksi
  ixi1=xpn/dpn
  ixi2=ypn/dpn
  ixi3=zpn/dpn

  ! Unit vectors of the face normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are

  ! Angle between vectorsa n and i_xi - we need cosine
  costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

  ! Relaxation factor for higher-order cell face gradient
  ! In general, nrelax can be any signed integer from some 
  ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
  !costn = costheta**nrelax

  costn = 1.0_dp

  if(nrelax == 1) then
    ! Minimal correction: nrelax = +1 :
    costn = costheta
  elseif(nrelax == 0) then
    ! Orthogonal correction: nrelax =  0 : 
    costn = 1.0_dp
  elseif(nrelax == -1) then
    ! Over-relaxed approach: nrelax = -1 :
    costn = 1.0_dp/costheta  
  endif

  ! dpp_j * sf
  vole=xpn*arx+ypn*ary+zpn*arz


  ! Interpolate gradients defined at CV centers to faces
  dfixi = dFidxi(1,ijp)*fxp+dFidxi(1,ijn)*fxn
  dfiyi = dFidxi(2,ijp)*fxp+dFidxi(2,ijn)*fxn
  dfizi = dFidxi(3,ijp)*fxp+dFidxi(3,ijn)*fxn


  !-- Skewness correction -->
  if (adjustl(trim(approach)) == 'skewness') then

    ! Overrelaxed correction vector d2, where s=dpn+d2
    d1x = costn
    d1y = costn
    d1z = costn

    d2x = xpn*costn
    d2y = ypn*costn
    d2z = zpn*costn

    !.....du/dx_i interpolated at cell face:
    dfixii = dfixi*d1x + arx/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
    dfiyii = dfiyi*d1y + ary/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
    dfizii = dfizi*d1z + arz/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z )

  ! |-- Intersection point offset and skewness correction -->


  elseif (adjustl(trim(approach)) == 'offset') then

    ! Find points P' and Pj'
    xpp=xf-(xf-xc(ijp))*nxx
    ypp=yf-(yf-yc(ijp))*nyy 
    zpp=zf-(zf-zc(ijp))*nzz

    xep=xf-(xf-xc(ijn))*nxx 
    yep=yf-(yf-yc(ijn))*nyy 
    zep=zf-(zf-zc(ijn))*nzz     

    xpnp = xep-xpp 
    ypnp = yep-ypp 
    zpnp = zep-zpp

    volep = arx*xpnp+ary*ypnp+arz*zpnp

   ! Overrelaxed correction vector d2, where S=dpn+d2
    d1x = costn
    d1y = costn
    d1z = costn
    
    xpnp = xpnp*costn
    ypnp = ypnp*costn
    zpnp = zpnp*costn

    ! The cell face interpolated gradient (d phi / dx_i)_j:
    ! Nonorthogonal corrections:          ___
    ! nablaFIxdnnp =>> dot_product(dFidxi,dNN')
    ! And:                                ___
    ! nablaFIxdnnp =>> dot_product(dFidxi,dPP')
    nablaFIxdnnp = dFidxi(1,ijn)*(xep-xc(ijn))+dFidxi(2,ijn)*(yep-yc(ijn))+dFidxi(3,ijn)*(zep-zc(ijn))
    nablaFIxdppp = dFidxi(1,ijp)*(xpp-xc(ijp))+dFidxi(2,ijp)*(ypp-yc(ijp))+dFidxi(3,ijp)*(zpp-zc(ijp))

    dfixii = dfixi*d1x + arx/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
    dfiyii = dfiyi*d1y + ary/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
    dfizii = dfizi*d1z + arz/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
 
  !-- Uncorrected -->
  elseif (adjustl(trim(approach)) == 'uncorrected') then
    
    dfixii = dfixi
    dfiyii = dfiyi
    dfizii = dfizi

  endif


end subroutine


!******************************************************************************
!
subroutine sngrad_vector_field(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
                               u,v,w, dudxi,dvdxi,dwdxi, nrelax, approach, &
                               duxi, duyi, duzi, dvxi, dvyi, dvzi, dwxi, dwyi, dwzi, &
                               duxii, dvxii, dwxii, duyii, dvyii, dwyii, duzii, dvzii, dwzii)
!
!******************************************************************************
!
!  Surface normal gradient with non-orthogonal correction done in two
!  possible ways - either by skewness correction of intersection point
!  offset.
!
!  Check out reference paper: 
!    Mirkov, Rasuo, Kenjeres, JCP, Vol. 287, 2015.
!
! Note: same as above, just for vector field.
!
!******************************************************************************
!
  implicit none
!
!******************************************************************************
! 
  integer, intent(in) :: ijp, ijn
  integer, intent(in) :: nrelax
  character(len=12) :: approach
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: lambda
  real(dp), dimension(numTotal), intent(in) :: u,v,w
  real(dp), dimension(3,numPCells), intent(in) :: dudxi,dvdxi,dwdxi
  real(dp), intent(out) ::  duxi,duyi,duzi,dvxi,dvyi,dvzi,dwxi,dwyi,dwzi
  real(dp), intent(out) ::  duxii,dvxii,dwxii,duyii,dvyii,dwyii,duzii,dvzii,dwzii

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              u, dudxi, nrelax, approach, duxi, duyi, duzi, duxii, duyii, duzii)

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              v, dvdxi, nrelax, approach, dvxi, dvyi, dvzi, dvxii, dvyii, dvzii)

  call sngrad(ijp, ijn, xf, yf, zf, arx, ary, arz, lambda, &
              w, dwdxi, nrelax, approach, dwxi, dwyi, dwzi, dwxii, dwyii, dwzii)

end subroutine


end module gradients