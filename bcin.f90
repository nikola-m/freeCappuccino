subroutine bcin

  use types
  use parameters
  use geometry
  use variables
  use title_mod, only: inlet_file
  use k_epsilon_std, only: te,ed ! mozda ovo drugacije, sta cemo kad je neki drugi model
  use temperature, only: t

  implicit none

  integer :: i,ini,ino
  real(dp) :: uav,outare
  real(dp) :: flowen, flowte, flowed
  real(dp) :: are


  flomas = 0.0_dp 
  flomom = 0.0_dp
  flowen = 0.0_dp
  flowte = 0.0_dp
  flowed = 0.0_dp

  open(unit=7,file=inlet_file)
  rewind 7

  ! Loop over inlet boundaries
  do i = 1,ninl
    ini = iInletStart + i

    read(7,*) u(ini),v(ini),w(ini),p(ini),te(ini),ed(ini),t(ini)

    vis(ini) = viscos

    fmi(i) = den(ini)*(Xni(i)*u(ini)+Yni(i)*v(ini)+Zni(i)*w(ini))

    ! Face normal vector is faced otwards, while velocity vector at inlet
    ! is faced inwards. That means their scalar product will be negative,
    ! so minus signs there is to turn net mass influx - flomas, into positive value.
    flomas = flomas - fmi(i) 

    flomom = flomom + fmi(i)*sqrt(u(ini)**2+v(ini)**2+w(ini)**2)

    flowen = flowen + abs(fmi(i)*t(ini))

    flowte = flowte + abs(fmi(i)*te(ini))

    flowed = flowed + abs(fmi(i)*ed(ini))

  enddo

  close(7)

  ! Loop over outlet boundaries

  ! Outlet area
  outare = 0.0_dp
  do i = 1,nout
    outare = outare + sqrt(Xno(i)**2+Yno(i)**2+Zno(i)**2)
  enddo

  ! Average velocity at outlet boundary
  uav = flomas/(densit*outare+small)

  ! Mass flow trough outlet faces using Uav velocity
  do i = 1,nout
    ino = iOutletStart + i

    are = sqrt(Xno(i)**2+Yno(i)**2+Zno(i)**2)

    u(ino) = uav * Xno(i)/are
    v(ino) = uav * Yno(i)/are
    w(ino) = uav * Zno(i)/are

    !fmo(i)=den(ino)*uav*(Xno(i)+Yno(i)+Zno(i)) ! Originalno u caffa kodu.
    fmo(i)=den(ino)*(Xno(i)*u(ino)+Yno(i)*v(ino)+Zno(i)*w(ino))

  enddo

  ! No inflow into the domain: eg. natural convection case, etc.
  if(flomas.lt.small) flomas=1.0_dp
  if(flomom.lt.small) flomom=1.0_dp


  ! Loop over wall boundaries to calculate normal distance from cell center dnw.
  ! do iface = 1,nwal
  !   inw = ijw(iface)
  !   ijp = ijpw(iface)

  !   ! Face area 
  !   Are=sqrt(Xnw(iface)**2+Ynw(iface)**2+Znw(iface)**2)

  !   ! Face normals
  !   nxf=Xnw(Iface)/(Are+Small)
  !   nyf=Ynw(Iface)/(Are+Small)
  !   nzf=Znw(Iface)/(Are+Small)

  !   ! Distanca od ijp do ijb projektovana na normalu 
  !   dnw(i) = (xc(ijp)-xfw(iface))*nxf + (yc(ijp)-yfw(iface))*nyf + (zc(ijp)-zfw(iface))*nzf
  ! enddo



  ! Initialization of residual for all variables
  do i=1,nphi
    rnor(i) = 1.0_dp
    resor(i) = 0.0_dp
  enddo

  rnor(iu)    = 1.0_dp/(flomom+small)
  rnor(iv)    = rnor(iu)
  rnor(iw)    = rnor(iu)

  rnor(ip)    = 1.0_dp/(flomas+small)
  rnor(ite)   = 1.0_dp/(flowte+small)
  rnor(ied)   = 1.0_dp/(flowed+small)


   ! Correct turbulence at inlet for appropriate turbulence model
  if(lturb) call correct_turbulence_inlet()

end subroutine