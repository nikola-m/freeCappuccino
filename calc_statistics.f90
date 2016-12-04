!***********************************************************************
!
subroutine calc_statistics
!
!***********************************************************************
!
  use types
  use parameters
  use geometry, only: numCells
  use variables
  use statistics

  implicit none
!
!***********************************************************************
!
  integer :: inp
  real(dp) :: u_nsample,v_nsample,w_nsample!,te_nsample,con_nsample


  n_sample = n_sample+1

  do inp = 1,numCells

    ! Velocity field
    u_aver(inp) = u_aver(inp)+u(inp) 
    v_aver(inp) = v_aver(inp)+v(inp) 
    w_aver(inp) = w_aver(inp)+w(inp) 
    ! con_aver(inp) = con_aver(inp)+con(inp)

    ! Ensemble average over n samples
    u_nsample = u_aver(inp)/n_sample
    v_nsample = v_aver(inp)/n_sample
    w_nsample = w_aver(inp)/n_sample
    !      con_nsample = con_aver(inp)/n_sample

    ! Reynolds stress components
    uu_aver(inp) = uu_aver(inp)+(u(inp)-u_nsample)**2
    vv_aver(inp) = vv_aver(inp)+(v(inp)-v_nsample)**2
    ww_aver(inp) = ww_aver(inp)+(w(inp)-w_nsample)**2

    te_aver(inp) = (uu_aver(inp)+vv_aver(inp)+ww_aver(inp)) &
                      /(2*n_sample)

    uv_aver(inp) = uv_aver(inp)+ &
                 ((u(inp)-u_nsample)*(v(inp)-v_nsample))
    uw_aver(inp) = uw_aver(inp)+ &
                 ((u(inp)-u_nsample)*(w(inp)-w_nsample))
    vw_aver(inp) = vw_aver(inp)+ &
                 ((v(inp)-v_nsample)*(w(inp)-w_nsample))

!    concon_aver(inp) = concon_aver(inp)+ &
!                 (con(inp)-con_nsample)**2 

!    ucon_aver(inp) = ucon_aver(inp)+ &
!                 ((u(inp)-u_nsample)*(con(inp)-con_nsample))
!    vcon_aver(inp) = vcon_aver(inp)+ &
!                 ((v(inp)-v_nsample)*(con(inp)-con_nsample))
!    wcon_aver(inp) = wcon_aver(inp)+ &
!                 ((w(inp)-w_nsample)*(con(inp)-con_nsample))
  
  end do    

end subroutine
