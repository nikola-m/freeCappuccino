!***********************************************************************
!
subroutine dpcg(fi,ifi)
!
!***********************************************************************
!
!    This routine incorporates the diagonally preconditioned 
!    Conjugate Gradient solver for symmetric matrices in CSR sparse matrix format
!
!    Writen by nikola mirkov, 2016. nmirkov@vin.bg.ac.rs
!
!***********************************************************************
!
  use types 
  use parameters
  use geometry, only: numCells, numTotal, noc, ijl, ijr
  use sparse_matrix
  use title_mod

  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: ifi
  real(dp), dimension(numTotal), intent(inout) :: fi 

!
! Local variables
!
  integer :: i, k, ns, l
  real(dp), dimension(numCells) :: pk,zk
  real(dp) :: rsm, resmax, res0, resl
  real(dp) :: s0, sk, alf, bet, pkapk

! residual tolerance
  resmax = sor(ifi)

!
! Initalize working arrays
!
  pk = 0.0_dp
  zk = 0.0_dp
  res = 0.0_dp

!
! Calculate initial residual vector and the norm
!

! res(:) = su(:) - sum( a(ioffset(1:numCells):ioffset(2:numCells+1)-1) * fi(ja(ioffset(1:numCells):ioffset(2:numCells+1)-1)) )

  do i=1,numCells
    res(i) = su(i) 
    do k = ioffset(i),ioffset(i+1)-1
      res(i) = res(i) -  a(k) * fi(ja(k)) 
    enddo
  enddo

  do i=1,noc
    res(ijl(i)) = res(ijl(i)) - ar(i)*fi(ijr(i))
    res(ijr(i)) = res(ijr(i)) - al(i)*fi(ijl(i))
  end do

  ! L1-norm of residual
  res0=sum(abs(res))
    if( res0.lt.sor(ifi) ) return
!
! If ltest=true, print the norm 
!
  if(ltest) write(6,'(20x,a,1pe10.3)') 'res0 = ',res0

  s0=1.e20
!
! Start iterations
!
  ns=nsw(ifi)
  do l=1,ns
!
! Solve for zk(ijk) -- diagonal preconditioner
!
  do i=1,numCells
    zk(i) = res(i) / a( diag(i) )
  enddo
  
  ! Inner product
  sk = sum(res*zk) !..or  dot_product(res,zk)

!
! Calculate beta
!
  bet=sk/s0

!
! Calculate new search vector pk
!
  pk = zk + bet*pk

!
! Calculate scalar product (pk.a pk) and alpha (overwrite zk)
!
  do i=1,numCells
    zk(i) = 0.0_dp 
    do k = ioffset(i),ioffset(i+1)-1
      zk(i) = zk(i) + a(k) * pk( ja(k) ) 
    enddo
  enddo

  ! block-boundaries
  do i=1,noc
    zk(ijl(i)) = zk(ijl(i)) + ar(i)*pk(ijr(i))
    zk(ijr(i)) = zk(ijr(i)) + al(i)*pk(ijl(i))
  end do

  ! Inner product
  pkapk=sum(pk*zk)

  alf=sk/pkapk

  ! Update solution vector
  fi(1:numCells) = fi(1:numCells) + alf*pk

  ! Update residual vector
  res = res - alf*zk

  ! L^1-norm of residual
  resl = sum(abs(res))

  s0=sk

!
! Check convergence
!
  if(l.eq.1) resor(ifi) = res0
  rsm = resl/(resor(ifi)+small)
  if(ltest) write(6,'(19x,3a,i4,a,1pe10.3,a,1pe10.3)') ' fi=',chvar(ifi),' sweep = ',l,' resl = ',resl,' rsm = ',rsm
  if(rsm.lt.resmax) exit
!
! End of iteration loop
!
  end do

! Write linear solver report:
  write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') &
  'PCG(Jacobi):  Solving for ',trim(chvarSolver(IFI)), &
  ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L

end subroutine
