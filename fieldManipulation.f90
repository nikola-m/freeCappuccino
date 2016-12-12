    module fieldManipulation
!
!***********************************************************************
!
 
   public

    contains

    pure function volumeWeightedAverage(U) result(wAvgU)
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


!***********************************************************************
!
      pure function face_interpolated(u,dUdxi,inp,inn,xf,yf,zf,lambda) result(ue)
!
!***********************************************************************
!
!     Variable interpolated at cell face center with non-orthogonality 
!     correction.
!     This is broken down version of two sided interpolation with cell
!     values at two sides and corresponding gradients.
!
!***********************************************************************
!
      use types
      use parameters
      use geometry, only: numTotal,numCells,xc,yc,zc

      implicit none
!
!***********************************************************************
!

!.....Result
      real(dp) :: ue
!.....Arguments
      integer, intent(in) :: inp, inn
      real(dp), dimension(numTotal), intent(in) :: u
      real(dp), dimension(3,numCells), intent(in) :: dUdxi
      real(dp), intent(in) :: xf,yf,zf
      real(dp), intent(in) :: lambda
!.....Locals
      real(dp) :: xpn,ypn,zpn,xi,yi,zi,fxp,fxn


      fxn=lambda
      fxp=1.0d0-fxn

!.....Distance vector between cell centers
      xpn=xc(inn)-xc(inp)
      ypn=yc(inn)-yc(inp)
      zpn=zc(inn)-zc(inp)

!.....Coordinates of cell-face center - j

!.....Coordinates of intersection point - j'
      Xi=Xc(Inp)*Fxp+Xc(Inn)*Fxn
      Yi=Yc(Inp)*Fxp+Yc(Inn)*Fxn
      Zi=Zc(Inp)*Fxp+Zc(Inn)*Fxn

      Ue = 0.5*( (U(Inp)+U(Inn)) +                        &
                 (                                        &
                   (dUdxi(1,Inp)+dUdxi(1,Inn))*(Xf-Xi) +  &
                   (dUdxi(2,Inp)+dUdxi(2,Inn))*(Yf-Yi) +  &
                   (dUdxi(3,Inp)+dUdxi(3,Inn))*(Zf-Zi)    &
                 ) +                                      &
                 (                                        &
                    dUdxi(1,Inp)*Xpn*Fxp +                &
                    dUdxi(2,Inp)*Ypn*Fxp +                &
                    dUdxi(3,Inp)*Zpn*Fxp                  &
                  ) +                                     &
                  (                                       &
                    dUdxi(1,Inn)*Xpn*Fxn +                &
                    dUdxi(2,Inn)*Ypn*Fxn +                &
                    dUdxi(3,Inn)*Zpn*Fxn                  &
                  )                                       &
               )

      return
      end function

      double precision function face_value(inp,inn, xf, yf, zf, fi, gradfi)
!=======================================================================
!     Calculates face value using values of variables and their gradients
!     at multiple neighbours cell-centers,
!     in least-squares sense.
!=======================================================================
      use types
      use parameters
      use geometry, only: numTotal,numCells,xc,yc,zc

      implicit none

!     Input
      integer :: inp, inn
      real(dp) :: xf, yf, zf
      real(dp), dimension(numTotal) :: fi
      real(dp), dimension(3,numCells) :: gradfi

!     Locals
      real(dp) ::  phi_p, phi_n
      real(dp) :: xcp,ycp,zcp
      real(dp) :: xcn,ycn,zcn
      real(dp) :: gradfi_p_x,gradfi_p_y,gradfi_p_z
      real(dp) :: gradfi_n_x,gradfi_n_y,gradfi_n_z
      real(dp) :: nr
      real(dp) :: gradfidr



!.....Values at cell center's of neighbouring cells:
      phi_p = fi(inp)

      phi_n = fi(inn)

      xcp = xc(inp)
      ycp = yc(inp)
      zcp = zc(inp)

      xcn = xc(inn)
      ycn = yc(inn)
      zcn = zc(inn)

      gradfi_p_x = gradfi(1,inp)
      gradfi_p_y = gradfi(2,inp)
      gradfi_p_z = gradfi(3,inp)

      gradfi_n_x = gradfi(1,inn)
      gradfi_n_y = gradfi(2,inn)
      gradfi_n_z = gradfi(3,inn)

      nr = 0.50D0

!.....gradfixdr = (sum(gradphi_nb(i,:)*r_nb2f(i,:)), i=1,n)
!@      gradfidr=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) 
      gradfidr=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) &
              +gradfi_n_x*(xf-xcn)+gradfi_n_y*(yf-ycn)+gradfi_n_z*(zf-zcn)

!@      face_value = ( phi_p + gradfidr )
      face_value = nr*( phi_p + phi_n + gradfidr)

      end function



      double precision function slope_limited_face_value(inp, xf, yf, zf, fi, gradfi, fimax,fimin) result(face_value)
!=======================================================================
!     Calculates face value using values of variables and their gradients
!     at cell-center
!     Cell-centered gradient limited using slope limiter:
!     Wang modified Venkataktirshnan slope limiter
!     Ref.: Z. J. Wang. "A Fast Nested Multi-grid Viscous Flow Solver for Adaptive Cartesian/Quad Grids",
!     International Journal for Numerical Methods in Fluids. 33. 657–680. 2000.
!     The same slope limiter is used in Fluent.
!=======================================================================
      use types
      use parameters
      use sparse_matrix
      use geometry, only: numTotal,numCells,xc,yc,zc

      implicit none

!     Input
      integer :: inp
      real(dp) :: xf, yf, zf
      real(dp),dimension(numTotal) :: fi
      real(dp),dimension(3,numCells) :: gradfi
      real(dp) :: fimax,fimin! NOTE: fimax i fimin, su globalni max i min u polju.

!     Locals
      integer :: i
      real(dp) :: phi_p
      real(dp) :: gradfiXdr,slopelimit
      real(dp) :: deltam,deltap,epsi
      real(dp) :: phi_max,phi_min

!.....Values at cell center:
      phi_p = fi(inp)

!.....gradfixdr = (sum(gradphi_nb(i,:)*r_nb2f(i,:)), i=1,n)
      gradfiXdr=gradfi(1,inp)*(xf-xc(inp))+gradfi(2,inp)*(yf-yc(inp))+gradfi(3,inp)*(zf-zc(inp)) 

!.....Find unlimited value:
      face_value =  phi_p + gradfiXdr 

!:::::Define slope limiter:

!.....max and min values over current cell and neighbors
      phi_max = fi(ja( ioffset(inp) ))
      phi_min = fi(ja( ioffset(inp) ))

      do i=ioffset(inp)+1, ioffset(inp+1)-1
        phi_max = max( phi_max, fi(ja(i)) )
        phi_min = min( phi_max, fi(ja(i)) )      
      enddo
      ! phi_max = max( fi(ja(ioffset(inp))) : fi(ja(ioffset(inp+1)-1)) ) ! NOTE: Mora loop jer nisu contiguous u memoriji!!!!
      ! phi_min = min( fi(ja(ioffset(inp))) : fi(ja(ioffset(inp+1)-1)) )

      deltam = face_value - phi_p
      if (deltam .gt. 0.0d0) then
          deltap = phi_max-phi_p
      else
          deltap = phi_min-phi_p
      endif

!.....Original Venkatakrishnan K=[0,?], we take fixed K=0.05
!      epsi = (0.05*vol(inp))**3 
  
!.....Wang proposition for epsilon
      epsi = (0.05*( fimax-fimin ))**2 

      slopelimit = 1./(deltam+small) *((deltap+epsi)*deltam+2*deltam**2*deltap) &
                                     /(deltap**2+2*deltam**2+deltap*deltam+epsi+small)


      face_value =  phi_p + slopelimit*gradfiXdr 

      end function


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
