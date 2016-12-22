module k_epsilon_std
!
! Implementation of Standard k-epsilon two equation turbulence model.
!
  use types
  use parameters
  use geometry
  use variables

  implicit none

  ! Turbulence model constants
  real(dp), parameter :: CMU = 0.09_dp
  real(dp), parameter :: C1 = 1.44_dp
  real(dp), parameter :: C2 = 1.92_dp
  real(dp), parameter :: C3 = 1.44_dp
  real(dp), parameter :: sigma_k = 1.0_dp
  real(dp), parameter :: sigma_epsilon = 1.3_dp

  ! Derived constants
  real(dp), parameter :: CMU25 = sqrt(sqrt(cmu))
  real(dp), parameter :: CMU75 = cmu25**3

  ! Unique identifier - here or in modules_allocatable???
  ! integer, parameter :: ite=5
  ! integer, parameter :: ied=6

  real(dp), dimension(:), allocatable :: te
  real(dp), dimension(:), allocatable :: ed

  real(dp), dimension(:), allocatable :: teo
  real(dp), dimension(:), allocatable :: edo

  real(dp), dimension(:), allocatable :: teoo
  real(dp), dimension(:), allocatable :: edoo
  
  real(dp), dimension(:,:), allocatable :: dTEdxi
  real(dp), dimension(:,:), allocatable :: dEDdxi



  private 

  public :: te, ed,  teo, edo,  teoo, edoo
  public :: dTEdxi, dEDdxi

  public :: allocate_k_epsilon_std
  public :: correct_turbulence_k_epsilon_std
  public :: correct_turbulence_inlet_k_epsilon_std

contains


subroutine allocate_k_epsilon_std
  use parameters
  implicit none

      allocate( te(numTotal) )
      allocate( ed(numTotal) )

      allocate( teo(numTotal) ) 
      allocate(edo(numTotal) )

      if( bdf .and. btime.gt.0.99 ) then 
        allocate( teoo(numTotal) ) 
        allocate( edoo(numTotal) ) 
      endif 

      allocate( dTEdxi(3,numCells) ) 
      allocate( dEDdxi(3,numCells) ) 

end subroutine allocate_k_epsilon_std



subroutine correct_turbulence_k_epsilon_std()
!
! Main module routine to solve turbulence model equations and subsequently update effective viscosity
!
  use types
  use parameters
  use variables
  use gradients
  implicit none

  call calcsc(TE,dTEdxi,ite) ! Assemble and solve turbulence kinetic energy eq.
  call calcsc(ED,dEDdxi,ied) ! Assemble and solve dissipation rate of tke eq.
  call modify_mu_eff()

end subroutine correct_turbulence_k_epsilon_std



subroutine correct_turbulence_inlet_k_epsilon_std()
!
! Update effective viscosity at inlet
!
  implicit none

  call modify_mu_eff_inlet()

end subroutine correct_turbulence_inlet_k_epsilon_std



subroutine calcsc(Fi,dFidxi,ifi)
!
!
!
  use types
  use parameters
  use geometry
  use variables
  use sparse_matrix
  use gradients
  use temperature, only: t,utt,vtt,wtt
  implicit none

  integer, intent(in) :: ifi
  real(dp), dimension(numTotal) :: fi
  real(dp), dimension(3,numCells) :: dfidxi

!
! Local variables
!
  integer ::  i, k, inp, ijp, ijn, ijb, iface
  real(dp) :: gam, prtr, apotime, const, urfrs, urfms, &
              utp, vtp, wtp, utn, vtn, wtn, &
              genp, genn, sut, &
              uttbuoy, vttbuoy, wttbuoy
  real(dp) :: cap, can, suadd
  real(dp) :: magStrainSq
  real(dp) :: off_diagonal_terms
  real(dp) :: are,nxf,nyf,nzf,vnp,xtp,ytp,ztp,ut2,tau
  real(dp) :: viss


  ! Variable specific coefficients:
  gam=gds(ifi)
  if(ifi.eq.ite) prtr=1.0d0/sigma_k
  if(ifi.eq.ite) prtr=1.0d0/sigma_epsilon

! Calculate gradient: 
  call grad(fi,dfidxi)

! Initialize source arrays
  su(:)=0.0d0
  sp(:)=0.0d0

!
! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
  if(ifi.eq.ite) then

!=========================================================
! STANDARD PRODUCTION
!=========================================================
  do inp=1,numCells
        magStrainSq=magStrain(inp)*magStrain(inp)
        gen(inp)=abs(vis(inp)-viscos)*magStrainSq
  enddo


!
!=====================================
! VOLUME SOURCE TERMS 
!=====================================
  do inp=1,numCells

        genp=max(gen(inp),zero)
        genn=min(gen(inp),zero)

      ! Add production term to the rhs:
        su(inp)=genp*vol(inp)              
      ! Add destruction term to the lhs:
        sp(inp)=ed(inp)*den(inp)*vol(inp)/(te(inp)+small)
        sp(inp)=sp(inp)-genn*vol(inp)/(te(inp)+small)


      !
      !=====================================
      ! VOLUME SOURCE TERMS: buoyancy
      !=====================================
        if(lcal(ien).and.lbuoy) then
          
          ! When bouy activated we need the freshest utt,vtt,wtt - turbulent heat fluxes
          call calcheatflux 

          if(boussinesq.eq.1) then
             uttbuoy=-gravx*den(inp)*utt(inp)*vol(inp)*beta
             vttbuoy=-gravy*den(inp)*vtt(inp)*vol(inp)*beta
             wttbuoy=-gravz*den(inp)*wtt(inp)*vol(inp)*beta
          else ! if (boussinesq.eq.0)
             uttbuoy=-gravx*den(inp)*utt(inp)*vol(inp)/(t(inp)+273.)
             vttbuoy=-gravy*den(inp)*vtt(inp)*vol(inp)/(t(inp)+273.)
             wttbuoy=-gravz*den(inp)*wtt(inp)*vol(inp)/(t(inp)+273.)
          end if

          utp=max(uttbuoy,zero)
          vtp=max(vttbuoy,zero)
          wtp=max(wttbuoy,zero)
          utn=min(uttbuoy,zero)
          vtn=min(vttbuoy,zero)
          wtn=min(wttbuoy,zero)

          su(inp)=su(inp)+utp+vtp+wtp
          sp(inp)=sp(inp)-(utn+vtn+wtn)/(te(inp)+small)

        end if

        !
        !=====================================
        ! UNSTEADY TERM
        ! Three Level Implicit Time Integration Method:
        ! in case that BTIME=0. --> Implicit Euler
        !=====================================
        if(bdf) then
          apotime=den(inp)*vol(inp)/timestep
          sut=apotime*((1+btime)*teo(inp)-0.5*btime*teoo(inp))
          su(inp)=su(inp)+sut
          sp(inp)=sp(inp)+apotime*(1+0.5*btime)
        endif

      !.....End of TE volume source terms
  enddo

!****************************************
  elseif(ifi.eq.ied) then
!****************************************

!
!=====================================
! VOLUME SOURCE TERMS 
!=====================================
  do inp=1,numCells

        genp=max(gen(inp),zero)
        genn=min(gen(inp),zero)

        su(inp)=c1*genp*ed(inp)*vol(inp)/(te(inp)+small)
        sp(inp)=c2*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)
        sp(inp)=sp(inp)-c1*genn*vol(inp)/(te(inp)+small) 

      !
      !=====================================
      ! VOLUME SOURCE TERMS: Buoyancy
      !=====================================
        if(lcal(ien).and.lbuoy) then
          const=c3*den(inp)*ed(inp)*vol(inp)/(te(inp)+small)

          if(boussinesq.eq.1) then
             uttbuoy=-gravx*utt(inp)*const*beta
             vttbuoy=-gravy*vtt(inp)*const*beta
             wttbuoy=-gravz*wtt(inp)*const*beta
          else ! if(boussinesq.eq.0)
             uttbuoy=-gravx*utt(inp)*const/(t(inp)+273.15)
             vttbuoy=-gravy*vtt(inp)*const/(t(inp)+273.15)
             wttbuoy=-gravz*wtt(inp)*const/(t(inp)+273.15)
          end if

          utp=max(uttbuoy,zero)
          vtp=max(vttbuoy,zero)
          wtp=max(wttbuoy,zero)
          utn=min(uttbuoy,zero)
          vtn=min(vttbuoy,zero)
          wtn=min(wttbuoy,zero)

          su(inp)=su(inp)+utp+vtp+wtp
          sp(inp)=sp(inp)-(utn+vtn+wtn)/(ed(inp)+small)
        end if

        !
        !=====================================
        !.....UNSTEADY TERM
        !=====================================
        if(bdf) then
        !
        !    Three Level Implicit Time Integration Method:
        !    in case that BTIME=0. --> Implicit Euler
        !
          apotime=den(inp)*vol(inp)/timestep
          sut=apotime*((1+btime)*edo(inp)-0.5*btime*edoo(inp))
          su(inp)=su(inp)+sut
          sp(inp)=sp(inp)+apotime*(1+0.5*btime)
        endif

      !.....End of IED volume source terms
  enddo
!--------------------------------------
  end if

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
        iface = iOCFacesStart+i
        ijp=ijl(i)
        ijn=ijr(i)
        call facefluxsc(ijp, ijn, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), fmoc(i), foc(i), gam, &
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
    iface = iInletFacesStart+i
    ijp = owner(iface)
    ijb = iInletStart+i

    ! dFidxi(:,ijb)=dFidxi(:,ijp) ! (constant gradient bc)

    call boundary_facefluxsc(ijp, ijb, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), fmi(i), &
     Fi, dFidxi, prtr, cap, can, suadd)

    Sp(ijp) = Sp(ijp)-can

    Su(ijp) = Su(ijp)-can*Fi(ijb)
  end do

  ! Outlet faces
  do i=1,nout
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i

    ! dFidxi(:,ijb)=dFidxi(:,ijp) ! (constant gradient bc)

    call boundary_facefluxsc(ijp, ijb, xf(iface), yf(iface), zf(iface), arx(iface), ary(iface), arz(iface), fmo(i), &
     FI, dFidxi, prtr, cap, can, suadd)

    Sp(ijp) = Sp(ijp)-can

    Su(ijp) = Su(ijp)-can*Fi(ijb)
  end do


  ! Wall boundary conditions
  if (ifi .eq. ite) then

    do i=1,nwal
      iface = iWallFacesStart+i
      ijp=owner(iface)
      ijb=iWallStart+i

      su(ijp)=su(ijp)-gen(ijp)*vol(ijp) ! oduzmi produkciju iz wall ajdecent celije
      viss=viscos
      if(ypl(i).gt.ctrans) viss=visw(i)

      ! Face area 
      are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

      ! Face normals
      nxf = arx(iface)/are
      nyf = ary(iface)/are
      nzf = arz(iface)/are

      ! Magnitude of a cell center velocity projected on boundary face normal
      Vnp = U(ijp)*nxf+V(ijp)*nyf+W(ijp)*nzf

      ! Tangential velocity components 
      xtp = U(ijp)-Vnp*nxf
      ytp = V(ijp)-Vnp*nyf
      ztp = W(ijp)-Vnp*nzf

      ! Its magnitude
      Vtp = xtp*xtp+ytp*ytp+ztp*ztp

      ! Tangent direction - unit vector
      xtp = xtp/vtp
      ytp = ytp/vtp
      ztp = ztp/vtp

      ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
      Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

      Tau = viss*Ut2/dnw(i)

      ! TAU=VISS*((U(IJB)-U(IJP))*XTW(IW) &
      !          +(V(IJB)-V(IJP))*YTW(IW))/DN(IW)

      gen(ijp)=abs(tau)*cmu25*sqrt(max(zero,te(ijp)))/(dnw(i)*cappa)
      su(ijp)=su(ijp)+gen(ijp)*vol(ijp)
    end do

  else

    ! Wall boundaries approximated with wall functions
    ! for correct values of dissipation all coefficients have
    ! to be zero, su equal the dissipation, and ap = 1

    do i=1,nwal
      iface = iWallFacesStart+i
      ijp=owner(iface)
      ijb=iWallStart+i

      ed(ijp)=cmu75*(max(zero,te(ijp)))**1.5/(cappa*dnw(i))
      su(ijp)=ed(ijp)

      a( ioffset(ijp):ioffset(ijp+1)-1 ) = 0.0_dp
      a( diag(ijp) ) = 1.0_dp

    end do

  endif

  ! Modify coefficients for Crank-Nicolson
  if (cn) then

      a(:) = 0.5_dp*a(:) ! Doesn't affect the main diagonal because it's still zero.

      if(ifi.eq.ite) then
        do i = 1,numInnerFaces
            ijp = owner(i)
            ijn = neighbour(i)

            k = icell_jcell_csr_value_index(i)
            su(ijp) = su(ijp) - a(k)*teo(ijn)

            k = jcell_icell_csr_value_index(i)
            su(ijn) = su(ijn) - a(k)*teo(ijp)
        enddo
        do ijp=1,numCells
            apotime=den(ijp)*vol(ijp)/timestep
            off_diagonal_terms = sum( a( ioffset(ijp) : ioffset(ijp+1)-1 ) ) !- a(diag(ijp))
            su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*teo(ijp)
            sp(ijp) = sp(ijp)+apotime
        enddo
      else ! ifi.eq.ied
        do i = 1,numInnerFaces
            ijp = owner(i)
            ijn = neighbour(i)

            k = icell_jcell_csr_value_index(i)
            su(ijp) = su(ijp) - a(k)*edo(ijn)

            k = jcell_icell_csr_value_index(i)
            su(ijn) = su(ijn) - a(k)*edo(ijp)
        enddo
        do ijp=1,numCells
            apotime=den(ijp)*vol(ijp)/timestep
            off_diagonal_terms = sum( a( ioffset(ijp) : ioffset(ijp+1)-1 ) ) !- a(diag(ijp))
            su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*edo(ijp)
            sp(ijp) = sp(ijp)+apotime
        enddo
      endif
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
    iface = iSymmetryFacesStart+i
    ijp = owner(iface)
    ijb = iSymmetryStart+i
    fi(ijb)=fi(ijp)
  end do
  !
 ! Outlet faces
  do i=1,nout
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i
    fi(ijb)=fi(ijp)
  end do

! These field values cannot be negative
  if(ifi.eq.ite.or.ifi.eq.ied) then
    fi(:)=max(fi(:),small)
  endif

end subroutine calcsc

!***********************************************************************
!
subroutine facefluxsc(ijp, ijn, xf, yf, zf, arx, ary, arz, flmass, lambda, gam, FI, dFidxi, prtr, cap, can, suadd)
!
!***********************************************************************
!
  use types
  use parameters
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
  game = (viste*prtr+viscos)

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
    !@fii=fi(ijp)*fxp+fi(ijn)*fxn+dfixi*(xf-xi)+dfiyi*(yf-yi)+dfizi*(zf-zi)
    !fii = face_interpolated(fi,dfidxi,inp,idew,idns,idtb,fxp,fxn)

    ! Explicit second order convection 
    !@fcfie=fm*fii
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

subroutine modify_mu_eff()
!
! Update turbulent and effective viscosity.
!
  use types
  use parameters
  use geometry
  use variables
  implicit none

  integer :: i,inp
  integer :: iface, ijp,ijb
  real(dp) :: visold
  real(dp) :: nxf,nyf,nzf,are
  real(dp) :: Vnp,Vtp,xtp,ytp,ztp
  real(dp) :: Ut2,Tau,Utau,ck,viscw

!==============================================================================
! Loop trough cells 
!==============================================================================

  do inp=1,numCells
        ! Store old value
        visold=vis(inp)
        ! Update effective viscosity:
        ! \mu_{eff}=\mu+\mu_t; \mu_t = C_\mu * \frac{k^2}{\epsilon} for standard k-epsilon
        vis(inp)=viscos + den(inp)*cmu*te(inp)**2/(ed(inp)+small)
        ! Underelaxation
        vis(inp)=urf(ivis)*vis(inp)+(1.0_dp-urf(ivis))*visold
  enddo

!==============================================================================
! Loop trough boundary faces
!==============================================================================

  !----------------------------------------------------------------------------
  ! Wall boundaries - update Visw and Ypl
  !----------------------------------------------------------------------------
  do i=1,nwal  
    iface = iWallFacesStart+i
    ijp = owner(iface)
    ijb = iWallStart+i

    ! Face area 
    are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

    ! Face normals
    nxf = arx(iface)/are
    nyf = ary(iface)/are
    nzf = arz(iface)/are

    ! Magnitude of a cell center velocity projected on boundary face normal
    Vnp = U(ijp)*nxf+V(ijp)*nyf+W(ijp)*nzf

    ! Tangential velocity components 
    xtp = U(ijp)-Vnp*nxf
    ytp = V(ijp)-Vnp*nyf
    ztp = W(ijp)-Vnp*nzf

    ! Its magnitude
    Vtp = xtp*xtp+ytp*ytp+ztp*ztp

    ! Tangent direction
    xtp = xtp/vtp
    ytp = ytp/vtp
    ztp = ztp/vtp

    ! projektovanje razlike brzina na pravac tangencijalne brzine u cell centru ijp
    Ut2 = abs( (U(ijb)-U(ijp))*xtp + (V(ijb)-V(ijp))*ytp + (W(ijb)-W(ijp))*ztp )

    Tau = viscos*Ut2/dnw(i)
    Utau = sqrt( Tau / den(ijb) )
    ypl(i) = den(ijb)*Utau*dnw(i)/viscos

    ! ! Ima i ova varijanta u cisto turb. granicni sloj varijanti sa prvom celijom u log sloju
    ! ck = cmu25*sqrt(max(te(ijp),zero))
    ! ypl(i) = den(ijb)*ck*dnw(i)/viscos
    ! ! ...ovo je tehnicki receno ystar iliti y* a ne y+

    viscw = zero

    if(ypl(i) > ctrans) then
      viscw = ypl(i)*viscos*cappa/log(Elog*ypl(i))
    endif

    visw(i) = max(viscos,viscw)
    vis(ijb) = visw(i)
  enddo
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Symmetry
  !----------------------------------------------------------------------------
  do i=1,nsym
    iface = iSymmetryFacesStart+i
    ijp = owner(iface)
    ijb = iSymmetryStart+i

    Vis(ijb) = Vis(ijp)
  enddo
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Inlet
  !----------------------------------------------------------------------------
  do i=1,ninl
    iface = iInletFacesStart+i
    ijp = owner(iface)
    ijb = iInletStart+i

    Vis(ijb) = Vis(ijp)
  enddo
  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! Outlet
  !----------------------------------------------------------------------------
  do i=1,nout
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijb = iOutletStart+i

    Vis(ijb) = Vis(ijp)
  enddo
  !----------------------------------------------------------------------------

end subroutine modify_mu_eff


subroutine modify_mu_eff_inlet()
!
! Update turbulent and effective viscosity at inlet.
!
  use types
  use parameters
  use geometry,only:ninl,iInletStart
  use variables

  implicit none

  integer :: i,ini

    ! Loop over inlet boundaries
  do i = 1,ninl
    ini = iInletStart+i   

    ! Update effective viscosity:
    ! \mu_{eff}=\mu+\mu_t; \mu_t = C_\mu * \frac{k^2}{\epsilon} for standard k-epsilon

    vis(ini) = viscos+den(ini)*te(ini)**2*cmu/(ed(ini)+small)

  enddo

end subroutine modify_mu_eff_inlet


end module k_epsilon_std






!***********************************************************************
!
subroutine boundary_facefluxsc(ijp, ijn, xf, yf, zf, arx, ary, arz, flmass, FI, dFidxi, prtr, cap, can, suadd)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: xc,yc,zc,numCells,numTotal
  use variables, only: vis

  implicit none
!
!***********************************************************************
! 

  integer, intent(in) :: ijp, ijn
  real(dp), intent(in) :: xf,yf,zf
  real(dp), intent(in) :: arx, ary, arz
  real(dp), intent(in) :: flmass
  ! real(dp), intent(in) :: lambda
  ! real(dp), intent(in) :: gam 
  real(dp), dimension(numTotal), intent(in) :: Fi
  real(dp), dimension(3,numTotal), intent(in) :: dFidxi
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

  real(dp) :: d1x,d1y,d1z,d2x,d2y,d2z

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
  ! fxn=lambda 
  ! fxp=1.0d0-lambda
  fxn=1.0_dp
  fxp=0.0_dp

  ! Distance vector between cell centers
  ! xpn=xc(ijn)-xc(ijp)
  ! ypn=yc(ijn)-yc(ijp)
  ! zpn=zc(ijn)-zc(ijp)
  xpn=xf-xc(ijp)
  ypn=yf-yc(ijp)
  zpn=zf-zc(ijp)

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
  game = (viste*prtr+viscos)




  ! Coordinates of point j'
  ! xi=xc(ijp)*fxp+xc(ijn)*fxn
  ! yi=yc(ijp)*fxp+yc(ijn)*fxn
  ! zi=zc(ijp)*fxp+zc(ijn)*fxn
  xi = xf
  yi = yf
  zi = zf



 !  !-- Intersection point offset and skewness correction --

 !  ! Find points P' and Pj'
 !  xpp=xf-(xf-xc(ijp))*nxx; ypp=yf-(yf-yc(ijp))*nyy; zpp=zf-(zf-zc(ijp))*nzz
 !  xep=xf-(xf-xc(ijn))*nxx; yep=yf-(yf-yc(ijn))*nyy; zep=zf-(zf-zc(ijn))*nzz     

 !  xpnp = xep-xpp; ypnp = yep-ypp; zpnp = zep-zpp
 !  volep = arx*xpnp+ary*ypnp+arz*zpnp

 !  ! Overrelaxed correction vector d2, where S=dpn+d2
 !  d1x = costn
 !  d1y = costn
 !  d1z = costn
  
 !  xpnp = xpnp*costn
 !  ypnp = ypnp*costn
 !  zpnp = zpnp*costn

 !  ! Interpolate gradients defined at CV centers to faces
 !  dfixi = dFidxi(1,ijp)*fxp+dFidxi(1,ijn)*fxn
 !  dfiyi = dFidxi(2,ijp)*fxp+dFidxi(2,ijn)*fxn
 !  dfizi = dFidxi(3,ijp)*fxp+dFidxi(3,ijn)*fxn

 !  ! The cell face interpolated gradient (d phi / dx_i)_j:
 !  ! Nonorthogonal corrections:         ___
 !  ! nablaFIxdnnp =>> dot_product(dFidxi,dNN')
 !  ! And:                               ___
 !  ! nablaFIxdnnp =>> dot_product(dFidxi,dPP')
 !  nablaFIxdnnp = dFidxi(1,ijn)*(xep-xc(ijn))+dFidxi(2,ijn)*(yep-yc(ijn))+dFidxi(3,ijn)*(zep-zc(ijn))
 !  nablaFIxdppp = dFidxi(1,ijp)*(xpp-xc(ijp))+dFidxi(2,ijp)*(ypp-yc(ijp))+dFidxi(3,ijp)*(zpp-zc(ijp))

 !  dfixii = dfixi*d1x + arx/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
 !  dfiyii = dfiyi*d1y + ary/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 
 !  dfizii = dfizi*d1z + arz/volep*( fi(ijn)+nablaFIxdnnp-fi(ijp)-nablaFixdppp-dfixi*xpnp-dfiyi*ypnp-dfizi*zpnp ) 

 !  !-- Intersection point offset and skewness correction --



  !-- Skewness correction --

  ! Overrelaxed correction vector d2, where s=dpn+d2
  d1x = costn
  d1y = costn
  d1z = costn

  d2x = xpn*costn
  d2y = ypn*costn
  d2z = zpn*costn

  ! Interpolate gradients defined at CV centers to faces
  ! dfixi = dFidxi(1,ijp)*fxp+dFidxi(1,ijn)*fxn
  ! dfiyi = dFidxi(2,ijp)*fxp+dFidxi(2,ijn)*fxn
  ! dfizi = dFidxi(3,ijp)*fxp+dFidxi(3,ijn)*fxn
  dfixi = dFidxi(1,ijp)
  dfiyi = dFidxi(2,ijp)
  dfizi = dFidxi(3,ijp) !...because constant gradient

  !.....du/dx_i interpolated at cell face:
  dfixii = dfixi*d1x + arx/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
  dfiyii = dfiyi*d1y + ary/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 
  dfizii = dfizi*d1z + arz/vole*( fi(ijn)-fi(ijp)-dfixi*d2x-dfiyi*d2y-dfizi*d2z ) 

  !-- Skewness correction --
 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

  ! if(lcds) then
  !   !---------------------------------------------
  !   ! CENTRAL DIFFERENCING SCHEME (CDS) 
  !   !---------------------------------------------
  !   ! Interpolate variable FI defined at CV centers to face using corrected CDS:
  !   !   |________Ue'___________|_______________Ucorr_____________________|
  !   fii=fi(ijp)*fxp+fi(ijn)*fxn+dfixi*(xf-xi)+dfiyi*(yf-yi)+dfizi*(zf-zi)
  !   !fii = face_interpolated(fi,dfidxi,inp,idew,idns,idtb,fxp,fxn)

  !   ! Explicit second order convection 
  !   fcfie=fm*fii
  ! else
  !   !---------------------------------------------
  !   ! Darwish-Moukalled TVD schemes for unstructured grids, IJHMT, 2003. 
  !   !---------------------------------------------
  !   ! Find r's - the gradient ratio. This is universal for all schemes.
  !   ! If flow goes from P to E
  !   r1 = (2*dFidxi(1,ijp)*xpn + 2*dFidxi(2,ijp)*ypn + 2*dFidxi(3,ijp)*zpn)/(FI(ijn)-FI(ijp)) - 1.0d0  
  !   ! If flow goes from E to P
  !   r2 = (2*dFidxi(1,ijn)*xpn + 2*dFidxi(2,ijn)*ypn + 2*dFidxi(3,ijn)*zpn)/(FI(ijp)-FI(ijn)) - 1.0d0 
  !   ! Find Psi for [ MUSCL ] :
  !   psiw = max(0., min(2.*r1, 0.5*r1+0.5, 2.))
  !   psie = max(0., min(2.*r2, 0.5*r2+0.5, 2.))
  !   ! High order flux at cell face
  !   fcfie =  ce*(fi(ijn) + fxn*psie*(fi(ijp)-fi(ijn)))+ &
  !            cp*(fi(ijp) + fxp*psiw*(fi(ijn)-fi(ijp)))
  ! endif

  ! ! Explicit first order convection
  ! fcfii = ce*fi(ijn)+cp*fi(ijp)

  ! Deffered correction for convection = gama_blending*(high-low)
  !ffic = gam*(fcfie-fcfii)

  !-------------------------------------------------------
  ! Explicit part of fluxes
  !-------------------------------------------------------
  ! suadd = -ffic+fdfie-fdfii 
  suadd = fdfie-fdfii 
  !-------------------------------------------------------

end subroutine boundary_facefluxsc