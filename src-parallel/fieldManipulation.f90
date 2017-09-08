    module fieldManipulation
!
!***********************************************************************
!
 
   public

    contains

    function volumeWeightedAverage(U) result(wAvgU)
    use types
    use parameters
    use geometry, only:numTotal,numCells,vol

    implicit none
!
!***********************************************************************
!
!...Output
    real(dp) :: wAvgU

!...Input
    real(dp), dimension(numTotal), intent(in) :: U

!...Locals
    integer :: inp
    real(dp) :: sumvol
  
    sumvol = 0.0_dp
    wAvgU = 0.0_dp 

      do inp=1,numCells
          wAvgU = wAvgU + (Vol(inp)*U(inp))
          sumvol = sumvol + vol(inp)
      enddo
    
    call global_sum( wAvgU )
    call global_sum( sumvol )

    wAvgU = wAvgU / sumvol

    end function



!***********************************************************************
!
      subroutine add_random_noise_to_field(Phi,percent)
!
!***********************************************************************
!
!     Level is Max perturbation aplitude in percent (%) from given field value.
!   
!     Example usage:
!       call add_random_noise_to_field(U,10)
!
      use types
      use parameters
      use geometry, only: numCells,numTotal

      implicit none
!
!***********************************************************************
!
      real(dp), dimension(numTotal), intent(inout) :: phi
      integer, intent(in) :: percent
      
      integer :: inp
      real(dp) :: level
      real(dp) :: perturb

      level = dble(percent)

      do inp=1,numCells

            ! Random number based fluctuation of mean profile            
            CALL init_random_seed()
            CALL RANDOM_NUMBER(perturb)
            
            ! perturb is now between 0. and 1., we want it to be from 0 to 2*amplitude
            ! e.g. perturb = 0.9+perturb/5. when Max perturbation is +/- 10% of mean profile
            perturb = ( 1.0_dp - level/100.0_dp ) + perturb * (2*level/100.0_dp)

            Phi(INP) = perturb*Phi(inp)

      enddo

      end subroutine


! !***********************************************************************
! !
! pure function von_karman_lengthscale() result(lvk)
! !
! !***********************************************************************
! !
!   use types
!   use parameters
!   use indexes
!   use geometry
!   use variables
!   use gradients
  
!   implicit none
! !
! !***********************************************************************
! !
! !....Output
!  real(dp), dimension(numCells) :: lvk

! !....Input
!  ! (None)

! !....Locals
!  integer :: inp
!  real(dp) :: uscnd
!  real(dp) :: d2udx2,d2udy2,d2udz2
!  real(dp) :: d2vdx2,d2vdy2,d2vdz2
!  real(dp) :: d2wdx2,d2wdy2,d2wdz2
 
!   do inp=1,numCells
     
!       volr = 1./vol(inp)
      

! !.....derivatives in  x- direction: 
! !      
!       dudxe=dUdxi(1,inp+nj)*fx(inp)+dUdxi(1,inp)*(1.0d0-fx(inp))
!       dudxw=dUdxi(1,inp)*fx(inp-nj)+dUdxi(1,inp-nj)*(1.0d0-fx(inp-nj))
!       dudxn=dUdxi(1,inp+1)*fy(inp)+dUdxi(1,inp)*(1.0d0-fy(inp))
!       dudxs=dUdxi(1,inp)*fy(inp-1)+dUdxi(1,inp-1)*(1.0d0-fy(inp-1))
!       dudxt=dUdxi(1,inp+nij)*fz(inp)+dUdxi(1,inp)*(1.0d0-fz(inp))
!       dudxb=dUdxi(1,inp)*fz(inp-nij)+dUdxi(1,inp-nij)*(1.0d0-fz(inp-nij))
      
!       dvdxe=dVdxi(1,inp+nj)*fx(inp)+dVdxi(1,inp)*(1.0d0-fx(inp))
!       dvdxw=dVdxi(1,inp)*fx(inp-nj)+dVdxi(1,inp-nj)*(1.0d0-fx(inp-nj))
!       dvdxn=dVdxi(1,inp+1)*fy(inp)+dVdxi(1,inp)*(1.0d0-fy(inp))
!       dvdxs=dVdxi(1,inp)*fy(inp-1)+dVdxi(1,inp-1)*(1.0d0-fy(inp-1))
!       dvdxt=dVdxi(1,inp+nij)*fz(inp)+dVdxi(1,inp)*(1.0d0-fz(inp))
!       dvdxb=dVdxi(1,inp)*fz(inp-nij)+dVdxi(1,inp-nij)*(1.0d0-fz(inp-nij))

!       dwdxe=dWdxi(1,inp+nj)*fx(inp)+dWdxi(1,inp)*(1.0d0-fx(inp))
!       dwdxw=dWdxi(1,inp)*fx(inp-nj)+dWdxi(1,inp-nj)*(1.0d0-fx(inp-nj))
!       dwdxn=dWdxi(1,inp+1)*fy(inp)+dWdxi(1,inp)*(1.0d0-fy(inp))
!       dwdxs=dWdxi(1,inp)*fy(inp-1)+dWdxi(1,inp-1)*(1.0d0-fy(inp-1))
!       dwdxt=dWdxi(1,inp+nij)*fz(inp)+dWdxi(1,inp)*(1.0d0-fz(inp))
!       dwdxb=dWdxi(1,inp)*fz(inp-nij)+dWdxi(1,inp-nij)*(1.0d0-fz(inp-nij))
! !
! !.....derivatives in y- direction:    
!       dudye=dUdxi(2,inp+nj)*fx(inp)+dUdxi(2,inp)*(1.0d0-fx(inp))
!       dudyw=dUdxi(2,inp)*fx(inp-nj)+dUdxi(2,inp-nj)*(1.0d0-fx(inp-nj))
!       dudyn=dUdxi(2,inp+1)*fy(inp)+dUdxi(2,inp)*(1.0d0-fy(inp))
!       dudys=dUdxi(2,inp)*fy(inp-1)+dUdxi(2,inp-1)*(1.0d0-fy(inp-1))
!       dudyt=dUdxi(2,inp+nij)*fz(inp)+dUdxi(2,inp)*(1.0d0-fz(inp))
!       dudyb=dUdxi(2,inp)*fz(inp-nij)+dUdxi(2,inp-nij)*(1.0d0-fz(inp-nij))

!       dvdye=dVdxi(2,inp+nj)*fx(inp)+dVdxi(2,inp)*(1.0d0-fx(inp))
!       dvdyw=dVdxi(2,inp)*fx(inp-nj)+dVdxi(2,inp-nj)*(1.0d0-fx(inp-nj))
!       dvdyn=dVdxi(2,inp+1)*fy(inp)+dVdxi(2,inp)*(1.0d0-fy(inp))
!       dvdys=dVdxi(2,inp)*fy(inp-1)+dVdxi(2,inp-1)*(1.0d0-fy(inp-1))
!       dvdyt=dVdxi(2,inp+nij)*fz(inp)+dVdxi(2,inp)*(1.0d0-fz(inp))
!       dvdyb=dVdxi(2,inp)*fz(inp-nij)+dVdxi(2,inp-nij)*(1.0d0-fz(inp-nij))

!       dwdye=dWdxi(2,inp+nj)*fx(inp)+dWdxi(2,inp)*(1.0d0-fx(inp))
!       dwdyw=dWdxi(2,inp)*fx(inp-nj)+dWdxi(2,inp-nj)*(1.0d0-fx(inp-nj))
!       dwdyn=dWdxi(2,inp+1)*fy(inp)+dWdxi(2,inp)*(1.0d0-fy(inp))
!       dwdys=dWdxi(2,inp)*fy(inp-1)+dWdxi(2,inp-1)*(1.0d0-fy(inp-1))
!       dwdyt=dWdxi(2,inp+nij)*fz(inp)+dWdxi(2,inp)*(1.0d0-fz(inp))
!       dwdyb=dWdxi(2,inp)*fz(inp-nij)+dWdxi(2,inp-nij)*(1.0d0-fz(inp-nij))
! !
! !.....derivatives in z- direction:      
!       dudze=dUdxi(3,inp+nj)*fx(inp)+dUdxi(3,inp)*(1.0d0-fx(inp))
!       dudzw=dUdxi(3,inp)*fx(inp-nj)+dUdxi(3,inp-nj)*(1.0d0-fx(inp-nj))
!       dudzn=dUdxi(3,inp+1)*fy(inp)+dUdxi(3,inp)*(1.0d0-fy(inp))
!       dudzs=dUdxi(3,inp)*fy(inp-1)+dUdxi(3,inp-1)*(1.0d0-fy(inp-1))
!       dudzt=dUdxi(3,inp+nij)*fz(inp)+dUdxi(3,inp)*(1.0d0-fz(inp))
!       dudzb=dUdxi(3,inp)*fz(inp-nij)+dUdxi(3,inp-nij)*(1.0d0-fz(inp-nij))

!       dvdze=dVdxi(3,inp+nj)*fx(inp)+dVdxi(3,inp)*(1.0d0-fx(inp))
!       dvdzw=dVdxi(3,inp)*fx(inp-nj)+dVdxi(3,inp-nj)*(1.0d0-fx(inp-nj))
!       dvdzn=dVdxi(3,inp+1)*fy(inp)+dVdxi(3,inp)*(1.0d0-fy(inp))
!       dvdzs=dVdxi(3,inp)*fy(inp-1)+dVdxi(3,inp-1)*(1.0d0-fy(inp-1))
!       dvdzt=dVdxi(3,inp+nij)*fz(inp)+dVdxi(3,inp)*(1.0d0-fz(inp))
!       dvdzb=dVdxi(3,inp)*fz(inp-nij)+dVdxi(3,inp-nij)*(1.0d0-fz(inp-nij))

!       dwdze=dWdxi(3,inp+nj)*fx(inp)+dWdxi(3,inp)*(1.0d0-fx(inp))
!       dwdzw=dWdxi(3,inp)*fx(inp-nj)+dWdxi(3,inp-nj)*(1.0d0-fx(inp-nj))
!       dwdzn=dWdxi(3,inp+1)*fy(inp)+dWdxi(3,inp)*(1.0d0-fy(inp))
!       dwdzs=dWdxi(3,inp)*fy(inp-1)+dWdxi(3,inp-1)*(1.0d0-fy(inp-1))
!       dwdzt=dWdxi(3,inp+nij)*fz(inp)+dWdxi(3,inp)*(1.0d0-fz(inp))
!       dwdzb=dWdxi(3,inp)*fz(inp-nij)+dWdxi(3,inp-nij)*(1.0d0-fz(inp-nij))
! !
! !.....second derivatives:
!       d2udx2 = ((dudxe-dudxw)*ar1x(inp)+(dudxn-dudxs)*ar2x(inp)+ &
!                 (dudxt-dudxb)*ar3x(inp))*volr 
!       d2udy2 = ((dudye-dudyw)*ar1y(inp)+(dudyn-dudys)*ar2y(inp)+ &
!                 (dudyt-dudyb)*ar3y(inp))*volr
!       d2udz2 = ((dudze-dudzw)*ar1z(inp)+(dudzn-dudzs)*ar2z(inp)+ &
!                 (dudzt-dudzb)*ar3z(inp))*volr
! !---------------
!       d2vdx2 = ((dvdxe-dvdxw)*ar1x(inp)+(dvdxn-dvdxs)*ar2x(inp)+ &
!                 (dvdxt-dvdxb)*ar3x(inp))*volr
!       d2vdy2 = ((dvdye-dvdyw)*ar1y(inp)+(dvdyn-dvdys)*ar2y(inp)+ &
!                 (dvdyt-dvdyb)*ar3y(inp))*volr
!       d2vdz2 = ((dvdze-dvdzw)*ar1z(inp)+(dvdzn-dvdzs)*ar2z(inp)+ &
!                 (dvdzt-dvdzb)*ar3z(inp))*volr
! !---------------
!       d2wdx2 = ((dwdxe-dwdxw)*ar1x(inp)+(dwdxn-dwdxs)*ar2x(inp)+ &
!                 (dwdxt-dwdxb)*ar3x(inp))*volr
!       d2wdy2 = ((dwdye-dwdyw)*ar1y(inp)+(dwdyn-dwdys)*ar2y(inp)+ &
!                 (dwdyt-dwdyb)*ar3y(inp))*volr
!       d2wdz2 = ((dwdze-dwdzw)*ar1z(inp)+(dwdzn-dwdzs)*ar2z(inp)+ &
!                 (dwdzt-dwdzb)*ar3z(inp))*volr
! !---------------


! !...2nd velocity derivative generalized to 3d using the magnitude of
! !   velocity laplacian
!     uscnd = sqrt((d2udx2+d2udy2+d2udz2)**2+ &
!                  (d2vdx2+d2vdy2+d2vdz2)**2+ &
!                  (d2wdx2+d2wdy2+d2wdz2)**2) 

! !.....von karman length scale
!     lvk(inp) = cappa*magStrain(inp)/uscnd
                 
!   end do

! end function

    
end module fieldManipulation
