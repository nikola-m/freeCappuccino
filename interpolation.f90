module interpolation


! private 

! public :: face_value

contains


! function face_value(ijp,ijn,xf,yf,zf,lambda,phi,dPhidxi,phimax,phimin) result(ue)

! use types
! use parameters

! implicit none

!   real(dp), dimension(numTotal), intent(in) :: phi
!   real(dp), dimension(3,numCells), intent(inout) :: dPhidxi

  
!   if () then 
!     ue = face_interpolated(ijp,ijn,xf,yf,zf,lambda,u,dUdxi)
!     ve = face_interpolated(ijp,ijn,xf,yf,zf,lambda,v,dVdxi)
!     we = face_interpolated(ijp,ijn,xf,yf,zf,lambda,w,dWdxi)
!   elseif () then 
!     ue = face_value_central(ijp, ijn, xf, yf, zf, u, dUdxi)
!     ve = face_value_central(ijp, ijn, xf, yf, zf, v, dVdxi)
!     we = face_value_central(ijp, ijn, xf, yf, zf, w, dWdxi)
!   elseif () then 
!     ue = face_value_2nd_upwind(ijp, xf, yf, zf, u, dUdxi)
!     ve = face_value_2nd_upwind(ijp, xf, yf, zf, v, dVdxi)
!     we = face_value_2nd_upwind(ijp, xf, yf, zf, w, dWdxi)
!   elseif () then
!     ue = face_value_2nd_upwind_slope_limited(ijp, xf, yf, zf, u, dUdxi, umin, umax)
!     ve = face_value_2nd_upwind_slope_limited(ijp, xf, yf, zf, v, dVdxi, vmin, vmax)
!     we = face_value_2nd_upwind_slope_limited(ijp, xf, yf, zf, w, dWdxi, wmin, wmax)
!   elseif () then
!     ue = face_value_muscl(ijp, ijn, xf, yf, zf, u, dUdxi)
!     ve = face_value_muscl(ijp, ijn, xf, yf, zf, v, dVdxi)
!     we = face_value_muscl(ijp, ijn, xf, yf, zf, w, dWdxi)
!   else
!     write(*,*) 'Interpolation scheme not chosen!'
!     stop
!   endif 

! end function




!***********************************************************************
!
      function face_value_central(inp,inn, xf, yf, zf, fi, gradfi) result(face_value)
!
!***********************************************************************
!
!     Calculates face value using values of variables and their gradients
!     at neighbours cell-centers.
!
!***********************************************************************
      use types
      use parameters
      use geometry, only: numTotal,numCells,xc,yc,zc

      implicit none

!.....Result
      real(dp) :: face_value

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

      nr = 0.5_dp
       
      gradfidr=gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp) &
              +gradfi_n_x*(xf-xcn)+gradfi_n_y*(yf-ycn)+gradfi_n_z*(zf-zcn)

      face_value = nr*( phi_p + phi_n + gradfidr)

      end function


!***********************************************************************
!
      function face_value_2nd_upwind(inp, xf, yf, zf, fi, gradfi) result(face_value)
!
!***********************************************************************
!
!     Calculates face value using values of variables and their gradients
!     at neighbours cell-centers.
!     Corresponds to unlimited second order upwind scheme as 
!     used in ANSYS FLUENT.
!
!***********************************************************************
!
      use types
      use parameters
      use geometry, only: numTotal,numCells,xc,yc,zc

      implicit none

!.....Result
      real(dp) :: face_value

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


      gradfidr = gradfi_p_x*(xf-xcp)+gradfi_p_y*(yf-ycp)+gradfi_p_z*(zf-zcp)

      face_value = phi_p + gradfidr

      end function


!***********************************************************************
!
      function face_value_muscl(inp,inn, xf, yf, zf, fi, gradfi) result(face_value)
!
!***********************************************************************
!
!     Calculates face value using values of variables and their gradients
!     at neighbours cell-centers.
!     Corresponds to MUSCL scheme as used in ANSYS FLUENT.
!
!***********************************************************************
!
      use types
      use parameters
      use geometry, only: numTotal,numCells,xc,yc,zc

      implicit none

!.....Result
      real(dp) :: face_value

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

      face_value = theta*face_value_central + (1.0_dp-theta)*face_value_2nd_upwind
      
      end function

!***********************************************************************
!
      function face_value_2nd_upwind_slope_limited(inp, xf, yf, zf, fi, gradfi, fimax, fimin) result(face_value)
!
!***********************************************************************
!
!     Calculates face value using values of variables and their gradients
!     at cell-center
!     Cell-centered gradient limited using slope limiter:
!     Wang modified Venkataktirshnan slope limiter
!     Ref.: Z. J. Wang. "A Fast Nested Multi-grid Viscous Flow Solver for Adaptive Cartesian/Quad Grids",
!     International Journal for Numerical Methods in Fluids. 33. 657–680. 2000.
!     The same slope limiter is used in Fluent.
!***********************************************************************
!
      use types
      use parameters
      use sparse_matrix, only: ioffset,ja
      use geometry, only: numTotal,numCells,xc,yc,zc

      implicit none

!     Result
      real(dp) :: face_value

!     Input
      integer :: inp
      real(dp) :: xf, yf, zf
      real(dp),dimension(numTotal) :: fi
      real(dp),dimension(3,numCells) :: gradfi
      real(dp) :: fimax,fimin! NOTE: fimax i fimin, su globalni max i min u polju.

!     Locals
      integer :: k
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

      do k=ioffset(inp)+1, ioffset(inp+1)-1
        phi_max = max( phi_max, fi(ja(k)) )
        phi_min = min( phi_max, fi(ja(k)) )      
      enddo

      deltam = face_value - phi_p
      if (deltam .gt. 0.0d0) then
          deltap = phi_max-phi_p
      else
          deltap = phi_min-phi_p
      endif

!.....Original Venkatakrishnan K=[0,?], we take fixed K=0.05
     ! epsi = (0.05*vol(inp))**3 
  
!.....Wang proposition for epsilon
      epsi = (0.05*( fimax-fimin ))**2 

      slopelimit = 1./(deltam+small) *((deltap+epsi)*deltam+2*deltam**2*deltap) &
                                     /(deltap**2+2*deltam**2+deltap*deltam+epsi+small)


      face_value =  phi_p + slopelimit*gradfiXdr 

      end function

end module