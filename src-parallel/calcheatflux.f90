!***********************************************************************
!
subroutine calcheatflux
!***********************************************************************
!
! Turbulent heat flux calculations.
!
! Three approaches are enabled:
! 1) Simple Gradient Diffusion Hypothesis (lsgdh = .true.)
! 2) Generalized Gradient Diffusion Hypothesis (lggdh = .true.)
! 3) Algebraic Flux Model (lafm = .true.)
!
! We need to set phit (default phit = 0.2) in input file for GGDH and AFM
! For all three we need to set underrelaxation factor facflx.
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix
  use variables
  use gradients, only: grad

  implicit none 
!
!***********************************************************************
!
  integer  :: inp

  real(dp) :: dtdx, dtdy, dtdz                    ! Temperature gradient vector components

  real(dp) ::  dudx, dudy, dudz, &                ! Velocity gradient
               dvdx, dvdy, dvdz, & 
               dwdx, dwdy, dwdz

  real(dp) :: volr, vist
  real(dp) :: tedi                                ! k/eps

  real(dp) ::  ut1, ut2, ut3, ut4, ut5, ut6, ut7
  real(dp) ::  vt1, vt2, vt3, vt4, vt5, vt6, vt7
  real(dp) ::  wt1, wt2, wt3, wt4, wt5, wt6, wt7
                          
  real(dp) :: uttold, vttold, wttold                


  ! Temperature gradient
  call grad(T,dTdxi)

  ! Velocity gradients of tentative velocity fields 
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

  do inp = 1,numCells

        volr = 1.0_dp/vol(inp)
        vist = (vis(inp)-viscos)/densit
        tedi = te(inp)/(ed(inp)+small)

        uttold = utt(inp)
        vttold = vtt(inp)
        wttold = wtt(inp)

        dTdx = dTdxi(1,inp)
        dTdy = dTdxi(2,inp)
        dTdz = dTdxi(3,inp)

    !
    !========================================
    !---[SGDH approach: ]
    !========================================
    if(lsgdh) then

        utt(inp) = -(vist*prt1)*dtdx
        vtt(inp) = -(vist*prt1)*dtdy
        wtt(inp) = -(vist*prt1)*dtdz
    !
    !========================================
    !---[GGDH approach: ]
    !========================================
    else if(lggdh) then

        ut1 = uu(inp)*dtdx
        ut2 = uv(inp)*dtdy
        ut3 = uw(inp)*dtdz

        vt1 = uv(inp)*dtdx
        vt2 = vv(inp)*dtdy
        vt3 = vw(inp)*dtdz

        wt1 = uw(inp)*dtdx
        wt2 = vw(inp)*dtdy
        wt3 = ww(inp)*dtdz

        utt(inp) = -phit*tedi*(ut1+ut2+ut3)
        vtt(inp) = -phit*tedi*(vt1+vt2+vt3)
        wtt(inp) = -phit*tedi*(wt1+wt2+wt3)
    !
    !========================================
    !.....[Algebraic Flux Model:]
    !========================================
    else if(lafm) then

    !--------------------------------
    !      [velocity gradients: ]
    !--------------------------------

        dudx = dUdxi(1,inp)
        dudy = dUdxi(2,inp)
        dudz = dUdxi(3,inp)

        dvdx = dVdxi(1,inp)
        dvdy = dVdxi(2,inp)
        dvdz = dVdxi(3,inp)

        dwdx = dWdxi(1,inp)
        dwdy = dWdxi(2,inp)
        dwdz = dWdxi(3,inp)

    !----------------------------------------------------------------
    !
        ut1 = uu(inp)*dtdx
        ut2 = uv(inp)*dtdy
        ut3 = uw(inp)*dtdz
        ut4 = sksi*utt(inp)*dudx
        ut5 = sksi*vtt(inp)*dudy
        ut6 = sksi*wtt(inp)*dudz
        ut7 = gravx*eta*vart(inp)*beta
        if(.not.boussinesq) ut7 = gravx*eta*vart(inp)/(t(inp)+273.0d0)

        vt1 = uv(inp)*dtdx
        vt2 = vv(inp)*dtdy
        vt3 = vw(inp)*dtdz
        vt4 = sksi*utt(inp)*dvdx
        vt5 = sksi*vtt(inp)*dvdy
        vt6 = sksi*wtt(inp)*dvdz
        vt7 = gravy*eta*vart(inp)*beta
        if(.not.boussinesq) vt7 = gravy*eta*vart(inp)/(t(inp)+273.0d0)

        wt1 = uw(inp)*dtdx
        wt2 = vw(inp)*dtdy
        wt3 = ww(inp)*dtdz
        wt4 = sksi*utt(inp)*dwdx
        wt5 = sksi*vtt(inp)*dwdy
        wt6 = sksi*wtt(inp)*dwdz
        wt7 = gravz*eta*vart(inp)*beta
        if(.not.boussinesq) wt7 = gravz*eta*vart(inp)/(t(inp)+273.0d0)

        utt(inp) = -phit*tedi*(ut1+ut2+ut3+ut4+ut5+ut6+ut7)
        vtt(inp) = -phit*tedi*(vt1+vt2+vt3+vt4+vt5+vt6+vt7)
        wtt(inp) = -phit*tedi*(wt1+wt2+wt3+wt4+wt5+wt6+wt7)

    !----------------------------------------------------------------
    !

    end if !! [conditional for SGDH,GGDH,AFM]


    ! Under-relaxation:
    utt(inp) = facflx*utt(inp)+(1.0_dp-facflx)*uttold
    vtt(inp) = facflx*vtt(inp)+(1.0_dp-facflx)*vttold
    wtt(inp) = facflx*wtt(inp)+(1.0_dp-facflx)*wttold

  end do ! cell-loop

end subroutine
