module interpolation

public

contains

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

      double precision function face_value_central(inp,inn, xf, yf, zf, fi, gradfi)
!=======================================================================
!     Calculates face value using values of variables and their gradients
!     at neighbours cell-centers.
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
       
      gradfidr=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) &
              +gradfi_n_x*(xf-xcn)+gradfi_n_y*(yf-ycn)+gradfi_n_z*(zf-zcn)

      face_value_central = nr*( phi_p + phi_n + gradfidr)

      end function

      double precision function face_value_2nd_upwind(inp, xf, yf, zf, fi, gradfi)
!=======================================================================
!     Calculates face value using values of variables and their gradients
!     at neighbours cell-centers.
!     Corresponds to unlimited second order upwind scheme as 
!     used in ANSYS FLUENT.
!=======================================================================
      use types
      use parameters
      use geometry, only: numTotal,numCells,xc,yc,zc

      implicit none

!     Input
      integer :: inp
      real(dp) :: xf, yf, zf
      real(dp), dimension(numTotal) :: fi
      real(dp), dimension(3,numCells) :: gradfi

!     Locals
      real(dp) ::  phi_p
      real(dp) :: xcp,ycp,zcp
      real(dp) :: gradfi_p_x,gradfi_p_y,gradfi_p_z
      real(dp) :: gradfidr



!.....Values at cell center's of neighbouring cells:
      phi_p = fi(inp)

      xcp = xc(inp)
      ycp = yc(inp)
      zcp = zc(inp)

      gradfi_p_x = gradfi(1,inp)
      gradfi_p_y = gradfi(2,inp)
      gradfi_p_z = gradfi(3,inp)


!.....gradfixdr = (sum(gradphi_nb(i,:)*r_nb2f(i,:)), i=1,n)
      gradfidr = gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp)

      face_value_2nd_upwind = phi_p + gradfidr

      end function


      double precision function face_value_muscl(inp,inn, xf, yf, zf, fi, gradfi)
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
      real(dp) :: gradfidr_2nd_upwind,gradfidr_central,face_value_2nd_upwind,face_value_central
      real(dp) :: theta

      ! theta = 1/8
      theta = 0.125_dp

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


!.....gradfixdr = (sum(gradphi_nb(i,:)*r_nb2f(i,:)), i=1,n)
      gradfidr_2nd_upwind=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) 
      gradfidr_central=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) &
                      +gradfi_n_x*(xf-xcn)+gradfi_n_y*(yf-ycn)+gradfi_n_z*(zf-zcn)

      face_value_2nd_upwind = ( phi_p + gradfidr_2nd_upwind )
      face_value_central = 0.5_dp*( phi_p + phi_n + gradfidr_central)

      face_value_muscl = theta*face_value_central + (1.0_dp-theta)*face_value_2nd_upwind
      
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



end module