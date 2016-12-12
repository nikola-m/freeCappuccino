subroutine find_strain_rate
!
!***********************************************************************
!
  use types
  use parameters  
  use geometry, only:numCells 
  use variables
  use gradients

  implicit none
!
!***********************************************************************
!
  integer :: inp
  real(dp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz 
  real(dp) :: s11,s12,s13,s21,s22,s23,s31,s32,s33,w12,w13,w23


  !.....Velocity gradients: 
  call grad(U,dUdxi)
  call grad(V,dVdxi)
  call grad(W,dWdxi)

  do inp=1,numCells

    dudx = dudxi(1,inp)
    dudy = dudxi(2,inp)
    dudz = dudxi(3,inp)

    dvdx = dvdxi(1,inp)
    dvdy = dvdxi(2,inp)
    dvdz = dvdxi(3,inp)

    dwdx = dwdxi(1,inp)
    dwdy = dwdxi(2,inp)
    dwdz = dwdxi(3,inp)

    ! Find strain rate tensor
    ! [s_ij]: |s_ij|=sqrt[2s_ij s_ij]
    s11=dudx
    s12=0.5*(dudy+dvdx)
    s13=0.5*(dudz+dwdx)
    s22=dvdy
    s23=0.5*(dvdz+dwdy) 
    s33=dwdz
    s21=s12  ! symmetry
    s31=s13
    s32=s23

    ! Find antisymmetric part of velocity gradient tensor
    ! [om_ij]: |om_ij|=sqrt[2 om_ij om_ij] 
    w12=0.5*(dudy - dvdx)
    w13=0.5*(dudz - dwdx)
    w23=0.5*(dvdz - dwdy)


    ! Find strain rate s = sqrt (2*sij*sij)
    magStrain(inp)=dsqrt(2*(s11**2+s22**2+s33**2 + 2*(s12**2+s13**2+s23**2)))

    ! Find scalar invariant of antisymmetric part of velocity gradient tensor
    vorticity(inp) = dsqrt(w12**2 + w23**2 + w13**2)

  enddo

end subroutine