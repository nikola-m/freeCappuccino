module temperature
!
! Implementation of sclar transport equation for temperature.
!
  use types
  use parameters
  use indexes
  use variables

  implicit none

  ! Constants
  real(dp), parameter :: sigt = 0.9_dp
  real(dp) :: pranl !(= 0.7_dp for air, 7.0_dp for water, read it from input file.)

 

  ! Unique identifier - here or in modules_allocatable???
  ! integer, parameter :: ien=7

  ! Temperature
  real(dp), dimension(:), allocatable :: t
  real(dp), dimension(:), allocatable :: to
  real(dp), dimension(:), allocatable :: too
  real(dp), dimension(:,:), allocatable :: dTdxi

  ! Temperature
  real(dp), dimension(:), allocatable :: vart
  real(dp), dimension(:), allocatable :: varto
  real(dp), dimension(:), allocatable :: vartoo
  real(dp), dimension(:,:), allocatable :: dVartdxi

  ! Heat fluxes
  real(dp), dimension(:), allocatable :: utt,vtt,wtt



  private 

  public :: t, to, too, vart, varto, vartoo, utt, vtt, wtt, pranl, &
            dTdxi, dVartdxi, &
            calcsc, &
            allocate_temperature, deallocate_temperature

contains


subroutine allocate_temperature
  use parameters
  implicit none

    allocate( t(numTotal) )
    allocate( to(numTotal) )

    if( bdf .and. btime.gt.0.5 ) then 
      allocate( too(numTotal) ) 
    endif 
    
    ! Gradient
    allocate( dTdxi(3,numCells) )

    ! Turbulent heat fluxes
    allocate( utt(numCells) )
    allocate( vtt(numCells) )
    allocate( wtt(numCells) )

end subroutine allocate_temperature



subroutine deallocate_temperature
  implicit none
    deallocate( t ) 
    deallocate( to )
    if (allocated(too)) deallocate( too )
    deallocate( dTdxi )
    deallocate( utt )
    deallocate( vtt )
    deallocate( wtt )
end subroutine deallocate_temperature



subroutine calcsc(Fi,dFidxi,ifi)
!
! Ansamble and solve transport eq. for temerature scalar.
!
  use types
  use parameters
  use indexes
  use geometry
  use variables
  use sparse_matrix
  use gradients
  implicit none

  integer, intent(in) :: ifi
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: dfidxi

!
! Local variables
!
  integer ::    i, k, inp, ijp, ijn, ijb
  real(dp) :: gam, prtr, apotime, urfrs, urfms
  real(dp) :: cap, can, suadd
  real(dp) :: off_diagonal_terms
  real(dp) :: coef,dcoef,sut


  ! Variable specific coefficients:
  gam=gds(ifi)
  prtr=1.0d0/sigt

! Calculate gradient: 
  call grad(fi,dfidxi)

! Initialize source arrays
  su(:)=0.0d0
  sp(:)=0.0d0

!
!=====================================
! VOLUME SOURCE TERMS 
!=====================================
  do inp=1,numCells

    !
    !=====================================
    ! UNSTEADY TERM
    ! Three Level Implicit Time Integration Method:
    ! in case that BTIME=0. --> Implicit Euler
    !=====================================
    if(bdf) then
      apotime=den(inp)*vol(inp)/timestep
      sut=apotime*((1+btime)*to(inp)-0.5*btime*too(inp))
      su(inp)=su(inp)+sut
      sp(inp)=sp(inp)+apotime*(1+0.5*btime)
    endif

    if(lturb.and.lbuoy) then
      ! When bouy activated we need the freshest utt,vtt,wtt - turbulent heat fluxes
      call calcheatflux 
      call Additional_algebraic_heatflux_terms
    end if

  enddo



! Calculate terms integrated over surfaces

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
        ijp = owner(i)
        ijn = neighbour(i)

        call facefluxsc(ijp, ijn, xf(i), yf(i), zf(i), arx(i), ary(i), arz(i), flmass(i), facint(i), gam, &
         fi, dFidxi, prtr, cap, can, suadd)

        ! > Off-diagonal elements:

        ! (icell,jcell) matrix element:
        k = icell_jcell_csr_value_index(i)
        a(k) = can

        ! (jcell,icell) matrix element:
        k = jcell_icell_csr_value_index(i)
        a(k) = cap

        ! > Elements on main diagonal:

        ! ! (icell,icell) main diagonal element
        ! k = diag(ijp)
        ! a(k) = a(k) - can
        sp(ijp) = sp(ijp) - can

        ! ! (jcell,jcell) main diagonal element
        ! k = diag(ijn)
        ! a(k) = a(k) - cap
        sp(ijn) = sp(ijn) - cap

        ! > Sources:

        su(ijp) = su(ijp) + suadd
        su(ijn) = su(ijn) - suadd 

  enddo


  ! Contribution from o- and c-grid cuts
  do i=1,noc
        ijp=ijl(i)
        ijn=ijr(i)
        call facefluxsc(ijp, ijn, xfoc(i), yfoc(i), zfoc(i), xnoc(i), ynoc(i), znoc(i), fmoc(i), foc(i), gam, &
         fi, dfidxi, prtr, al(i), ar(i), suadd)

        sp(ijp) = sp(ijp) - ar(i)
        sp(ijn) = sp(ijn) - al(i)

        su(ijp) = su(ijp) + suadd
        su(ijn) = su(ijn) - suadd

  end do

  !
  ! Boundary conditions
  !

  ! Inlet faces
  do i=1,ninl
    ijp = owner(iInletFacesStart+i)
    ijb = iInletStart+i

    dFidxi(:,ijb)=dFidxi(:,ijp) ! (constant gradient bc)

    call facefluxsc(ijp, ijb, xfi(i), yfi(i), zfi(i), xni(i), yni(i), zni(i), fmi(i), zero, one, &
     Fi, dFidxi, prtr, cap, can, suadd)

    Sp(ijp) = Sp(ijp)-can

    Su(ijp) = Su(ijp)-can*Fi(ijb)
  end do

  ! Outlet faces
  do i=1,nout
    ijp = owner(iOutletFacesStart+i)
    ijb = iOutletStart+i

    dFidxi(:,ijb)=dFidxi(:,ijp) ! (constant gradient bc)

    call facefluxsc(ijp, ijb, xfo(i), yfo(i), zfo(i), xno(i), yno(i), zno(i), fmo(i), zero, one, &
     FI, dFidxi, prtr, cap, can, suadd)

    Sp(ijp) = Sp(ijp)-can

    Su(ijp) = Su(ijp)-can*Fi(ijb)
  end do

  ! Wall boundary conditions

  ! Isothermal wall boundaries

  do i=1,nwali
    ijp=owner(iWallFacesStart+i)
    ijb=iWallStart+i
    dcoef = (viscos+(vis(ijp)-viscos)/sigt)/pranl ! Vrlo diskutabilno, proveriti!
    coef=dcoef*srdw(i)
    a(diag(ijp)) = a(diag(ijp)) + coef
    su(ijp) = su(ijp) + coef*t(ijb)
  end do

  ! Adiabatic wall boundaries

  do i=1,nwala
    ijp=owner(iWallFacesStart+nwali+i)
    ijb=iWallStart+nwali+i
    t(ijb)=t(ijp)
  end do

  ! Modify coefficients for Crank-Nicolson
  if (cn) then

      a(:) = 0.5_dp*a(:) ! Doesn't affect the main diagonal because it's still zero.

        do i = 1,numInnerFaces
            ijp = owner(i)
            ijn = neighbour(i)

            k = icell_jcell_csr_value_index(i)
            su(ijp) = su(ijp) - a(k)*to(ijn)

            k = jcell_icell_csr_value_index(i)
            su(ijn) = su(ijn) - a(k)*to(ijp)
        enddo
        do ijp=1,numCells
            apotime=den(ijp)*vol(ijp)/timestep
            off_diagonal_terms = sum( a( ioffset(ijp) : ioffset(ijp+1)-1 ) ) !- a(diag(ijp))
            su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*to(ijp)
            sp(ijp) = sp(ijp)+apotime
        enddo

  endif

  ! Underrelaxation factors
  urfrs=urfr(ifi)
  urfms=urfm(ifi)

  ! Main diagonal term assembly:
  do inp = 1,numCells
        ! Main diagonal term assembly:
        ! Sum all coefs in a row of a sparse matrix, but since we also included diagonal element 
        ! we substract it from the sum, to eliminate it from the sum.
        off_diagonal_terms  = sum( a(ioffset(inp) : ioffset(inp+1)-1) ) !- a(diag(ijp)) because = 0
        a(diag(inp)) = sp(inp) - off_diagonal_terms

        ! Underelaxation:
        a(diag(inp)) = a(diag(inp))*urfrs
        su(inp) = su(inp) + urfms*a(diag(inp))*fi(inp)

  enddo

  ! Solve linear system:
  call bicgstab(fi,ifi)

!
! Update symmetry and outlet boundaries
!
  ! Symmetry faces
  do i=1,nsym
    ijp = owner(iSymmetryFacesStart+i)
    ijb = iSymmetryStart+i
    fi(ijb)=fi(ijp)
  end do
  !
 ! Outlet faces
  do i=1,nout
    ijp = owner(iOutletFacesStart+i)
    ijb = iOutletStart+i
    fi(ijb)=fi(ijp)
  end do

! These field values cannot be negative
  fi(:)=max(fi(:),small)

end subroutine calcsc

!***********************************************************************
!
subroutine facefluxsc(ijp, ijn, xf, yf, zf, arx, ary, arz, flmass, lambda, gam, FI, dFidxi, prtr, cap, can, suadd)
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry, only: xc,yc,zc
  use variables, only: vis

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flmass
  real(dp), intent(in) :: lambda
  real(dp), intent(in) :: gam 
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numCells), intent(in) :: dFidxi
  real(dp), intent(in) :: prtr
  real(dp), intent(inout) :: cap, can, suadd


! Local variables
  real(dp) :: are
  real(dp) :: xpn,ypn,zpn
  real(dp) :: nxx,nyy,nzz,ixi1,ixi2,ixi3,dpn,costheta,costn
  real(dp) :: xi,yi,zi
  real(dp) :: Cp,Ce
  real(dp) :: fii,fm

  real(dp) :: fdfie,fdfii,fcfie,fcfii,ffic

  real(dp) :: d1x,d1y,d1z

  real(dp) :: de, vole, game, viste

  real(dp) :: fxp,fxn
  real(dp) :: xpp,ypp,zpp,xep,yep,zep,xpnp,ypnp,zpnp,volep
  real(dp) :: nablaFIxdnnp,nablaFIxdppp
  real(dp) :: dfixi,dfiyi,dfizi
  real(dp) :: dfixii,dfiyii,dfizii
  real(dp) :: r1,r2
  real(dp) :: psie,psiw
!----------------------------------------------------------------------

  dfixi = 0.0d0
  dfiyi = 0.0d0
  dfizi = 0.0d0

  ! > Geometry:

  ! Face interpolation factor
  fxn=lambda 
  fxp=1.0d0-lambda

  ! Distance vector between cell centers
  xpn=xc(ijn)-xc(ijp)
  ypn=yc(ijn)-yc(ijp)
  zpn=zc(ijn)-zc(ijp)

  ! Distance from P to neighbor N
  dpn=sqrt(xpn**2+ypn**2+zpn**2)     

  ! Components of the unit vector i_ksi
  ixi1=xpn/dpn
  ixi2=ypn/dpn
  ixi3=zpn/dpn

  ! cell face area
  are=sqrt(arx**2+ary**2+arz**2)

  ! Unit vectors of the normal
  nxx=arx/are
  nyy=ary/are
  nzz=arz/are

  ! Angle between vectorsa n and i_xi - we need cosine
  costheta=nxx*ixi1+nyy*ixi2+nzz*ixi3

  ! Relaxation factor for higher-order cell face gradient
  ! Minimal correction: nrelax = +1 :
  !costn = costheta
  ! Orthogonal correction: nrelax =  0 : 
  costn = 1.0d0
  ! Over-relaxed approach: nrelax = -1 :
  !costn = 1./costheta
  ! In general, nrelax can be any signed integer from some 
  ! reasonable interval [-nrelax,nrelax] (or maybe even real number): 
  !costn = costheta**nrelax

  ! dpp_j * sf
  vole=xpn*arx+ypn*ary+zpn*arz


  ! Cell face diffussion coefficint
  viste = (vis(ijp)-viscos)*fxp+(vis(ijn)-viscos)*fxn
  game = (viscos+viste*prtr)/pranl 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  ! Coordinates of point j'
  xi=xc(ijp)*fxp+xc(ijn)*fxn
  yi=yc(ijp)*fxp+yc(ijn)*fxn
  zi=zc(ijp)*fxp+zc(ijn)*fxn

  ! Find points P' and Pj'
  xpp=xf-(xf-xc(ijp))*nxx; ypp=yf-(yf-yc(ijp))*nyy; zpp=zf-(zf-zc(ijp))*nzz
  xep=xf-(xf-xc(ijn))*nxx; yep=yf-(yf-yc(ijn))*nyy; zep=zf-(zf-zc(ijn))*nzz     

  xpnp = xep-xpp; ypnp = yep-ypp; zpnp = zep-zpp
  volep = arx*xpnp+ary*ypnp+arz*zpnp

 ! Overrelaxed correction vector d2, where S=dpn+d2
  d1x = costn
  d1y = costn
  d1z = costn
  
  xpnp = xpnp*costn
  ypnp = ypnp*costn
  zpnp = zpnp*costn

  ! Interpolate gradients defined at CV centers to faces
  dfixi = dFidxi(1,ijp)*fxp+dFidxi(1,ijn)*fxn
  dfiyi = dFidxi(2,ijp)*fxp+dFidxi(2,ijn)*fxn
  dfizi = dFidxi(3,ijp)*fxp+dFidxi(3,ijn)*fxn

  ! The cell face interpolated gradient (d phi / dx_i)_j:
  ! Nonorthogonal corrections:         ___
  ! nablaFIxdnnp =>> dot_product(dFidxi,dNN')
  ! And:                               ___
  ! nablaFIxdnnp =>> dot_product(dFidxi,dPP')
  nablaFIxdnnp = dFidxi(1,ijn)*(xep-xc(ijn))+dFidxi(2,ijn)*(yep-yc(ijn))+dFidxi(3,ijn)*(zep-zc(ijn))
  nablaFIxdppp = dFidxi(1,ijp)*(xpp-xc(ijp))+dFidxi(2,ijp)*(ypp-yc(ijp))+dFidxi(3,ijp)*(zpp-zc(ijp))

  dfixii = dfixi*d1x + arx/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
  dfiyii = dfiyi*d1y + ary/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
  dfizii = dfizi*d1z + arz/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 


  ! Explicit diffusion
  fdfie = game*(dfixii*arx + dfiyii*ary + dfizii*arz)   
  ! Implicit diffussion 
  fdfii = game*are/dpn*(dfixi*xpn+dfiyi*ypn+dfizi*zpn)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  ! Difusion coefficient
  de = game*are/dpn

  ! Convection fluxes - uds
  fm = flmass
  ce = min(fm,zero) 
  cp = max(fm,zero)

  ! System matrix coefficients
  cap = -de - max(fm,zero)
  can = -de + min(fm,zero)

  if(lcds) then
    !---------------------------------------------
    ! CENTRAL DIFFERENCING SCHEME (CDS) 
    !---------------------------------------------
    ! Interpolate variable FI defined at CV centers to face using corrected CDS:
    !   |________Ue'___________|_______________Ucorr_____________________|
    fii=fi(ijp)*fxp+fi(ijn)*fxn+dfixi*(xf-xi)+dfiyi*(yf-yi)+dfizi*(zf-zi)
    !fii = face_interpolated(fi,dfidxi,inp,idew,idns,idtb,fxp,fxn)

    ! Explicit second order convection 
    fcfie=fm*fii
  else
    !---------------------------------------------
    ! Darwish-Moukalled TVD schemes for unstructured grids, IJHMT, 2003. 
    !---------------------------------------------
    ! Find r's - the gradient ratio. This is universal for all schemes.
    ! If flow goes from P to E
    r1 = (2*dFidxi(1,ijp)*xpn + 2*dFidxi(2,ijp)*ypn + 2*dFidxi(3,ijp)*zpn)/(FI(ijn)-FI(ijp)) - 1.0d0  
    ! If flow goes from E to P
    r2 = (2*dFidxi(1,ijn)*xpn + 2*dFidxi(2,ijn)*ypn + 2*dFidxi(3,ijn)*zpn)/(FI(ijp)-FI(ijn)) - 1.0d0 
    ! Find Psi for [ MUSCL ] :
    psiw = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
    psie = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
    ! High order flux at cell face
    fcfie =  ce*(fi(ijn) + fxn*psie*(fi(ijp)-fi(ijn)))+ &
             cp*(fi(ijp) + fxp*psiw*(fi(ijn)-fi(ijp)))
  endif

  ! Explicit first order convection
  fcfii = ce*fi(ijn)+cp*fi(ijp)
  ! Deffered correction for convection = gama_blending*(high-low)
  ffic = gam*(fcfie-fcfii)
  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  suadd = -ffic+fdfie-fdfii 
  !-------------------------------------------------------

end subroutine facefluxsc

end module temperature