!***********************************************************************
!
subroutine calcheatflux
!***********************************************************************
!
! Turbulent heat flux calculations  
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry
  use sparse_matrix
  use variables
  use gradients, only: grad
  use k_epsilon_std, only: te,ed
  use temperature, only: t,vart,dTdxi,utt,vtt,wtt

  implicit none 
!
!***********************************************************************
!
  integer :: inp
  real(dp) :: volr, vist, &
                tedi, &                                ! k/eps
                dtdx, dtdy, dtdz, &                    ! temperature gradient vector components
                dudx, dudy, dudz, &                    ! velocity gradient
                dvdx, dvdy, dvdz, & 
                dwdx, dwdy, dwdz , &
                ut1, ut2, ut3, ut4, ut5, ut6, ut7, &   ! needed for afm approach
                vt1, vt2, vt3, vt4, vt5, vt6, vt7, &
                wt1, wt2, wt3, wt4, wt5, wt6, wt7
                          
  real(dp) :: uttold, vttold, wttold                


  ! Temperature gradient
  call grad(T,dTdxi)
  
  ! if (lstsq) then 
  !   call grad_lsq_qr(T,dTdxi,2,d)
  ! else 
  !   call grad_gauss(T,dTdxi(1,:),dTdxi(2,:),dTdxi(3,:))
  ! endif

  ! Velocity gradients of tentative velocity fields 
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

  ! if (lstsq) then 
  !   call grad_lsq_qr(u,dUdxi,2,d)
  !   call grad_lsq_qr(v,dVdxi,2,d)
  !   call grad_lsq_qr(w,dWdxi,2,d)
  ! else
  !   call grad_gauss(u,dUdxi(1,:),dUdxi(2,:),dUdxi(3,:))
  !   call grad_gauss(v,dVdxi(1,:),dVdxi(2,:),dVdxi(3,:))
  !   call grad_gauss(w,dWdxi(1,:),dWdxi(2,:),dWdxi(3,:))
  ! endif

  do inp = 1,numCells

        volr=1./vol(inp)
        vist=(vis(inp)-viscos)/densit
        tedi=te(inp)/(ed(inp)+small)

        uttold=utt(inp)
        vttold=vtt(inp)
        wttold=wtt(inp)

        dTdx = dTdxi(1,inp)
        dTdy = dTdxi(2,inp)
        dTdz = dTdxi(3,inp)

    !
    !==========================================
    !---[SGDH approach: ]
    !==========================================
        if(lsgdh) then
        utt(inp)=-(vist*prt1)*dtdx
        vtt(inp)=-(vist*prt1)*dtdy
        wtt(inp)=-(vist*prt1)*dtdz
    !
    !==========================================
    !---[GGDH approach: ]
    !==========================================
        else if(lggdh) then
        ut1=uu(inp)*dtdx
        ut2=uv(inp)*dtdy
        ut3=uw(inp)*dtdz

        vt1=uv(inp)*dtdx
        vt2=vv(inp)*dtdy
        vt3=vw(inp)*dtdz

        wt1=uw(inp)*dtdx
        wt2=vw(inp)*dtdy
        wt3=ww(inp)*dtdz

        utt(inp)=-phit*tedi*(ut1+ut2+ut3)
        vtt(inp)=-phit*tedi*(vt1+vt2+vt3)
        wtt(inp)=-phit*tedi*(wt1+wt2+wt3)
    !
    !=======================================================
    !.....[AFM: include all]
    !=======================================================
          else if(lafm) then
    !--------------------------------
    !      [velocity gradients: ]
    !--------------------------------

        dudx=dUdxi(1,inp)
        dudy=dUdxi(2,inp)
        dudz=dUdxi(3,inp)

        dvdx=dVdxi(1,inp)
        dvdy=dVdxi(2,inp)
        dvdz=dVdxi(3,inp)

        dwdx=dWdxi(1,inp)
        dwdy=dWdxi(2,inp)
        dwdz=dWdxi(3,inp)

    !----------------------------------------------------------------
    !
        ut1=uu(inp)*dtdx
        ut2=uv(inp)*dtdy
        ut3=uw(inp)*dtdz
        ut4=sksi*utt(inp)*dudx
        ut5=sksi*vtt(inp)*dudy
        ut6=sksi*wtt(inp)*dudz
        ut7=gravx*eta*vart(inp)*beta
        if(boussinesq.eq.0)ut7=gravx*eta*vart(inp)/(t(inp)+273.0d0)

        vt1=uv(inp)*dtdx
        vt2=vv(inp)*dtdy
        vt3=vw(inp)*dtdz
        vt4=sksi*utt(inp)*dvdx
        vt5=sksi*vtt(inp)*dvdy
        vt6=sksi*wtt(inp)*dvdz
        vt7=gravy*eta*vart(inp)*beta
        if(boussinesq.eq.0)vt7=gravy*eta*vart(inp)/(t(inp)+273.0d0)

        wt1=uw(inp)*dtdx
        wt2=vw(inp)*dtdy
        wt3=ww(inp)*dtdz
        wt4=sksi*utt(inp)*dwdx
        wt5=sksi*vtt(inp)*dwdy
        wt6=sksi*wtt(inp)*dwdz
        wt7=gravz*eta*vart(inp)*beta
        if(boussinesq.eq.0)wt7=gravz*eta*vart(inp)/(t(inp)+273.0d0)

        utt(inp)=-phit*tedi*(ut1+ut2+ut3+ut4+ut5+ut6+ut7)
        vtt(inp)=-phit*tedi*(vt1+vt2+vt3+vt4+vt5+vt6+vt7)
        wtt(inp)=-phit*tedi*(wt1+wt2+wt3+wt4+wt5+wt6+wt7)

        end if !! [conditional for SGDH,GGDH,AFM]
    !-------------------------------------------------------
    !
        utt(inp)=facflx*utt(inp)+(1.0d0-facflx)*uttold
        vtt(inp)=facflx*vtt(inp)+(1.0d0-facflx)*vttold
        wtt(inp)=facflx*wtt(inp)+(1.0d0-facflx)*wttold

  end do ! cell-loop

end subroutine
