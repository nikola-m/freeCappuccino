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
  use geometry, only: numCells,numTotal,noc,ijl,ijr,iProcFacesStart,iProcStart
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
  integer :: i, k, ns, l, iOtherProc
  real(dp), dimension(numCells) :: zk
  real(dp), dimension(numCells+npro) :: pk
  real(dp) :: rsm, resmax, res0, resl, tol
  real(dp) :: s0, sk, alf, bet, pkapk

! residual tolerance
  resmax = sor(ifi)
  tol = 1e-13

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

  ! Residual contribution for cells at processor boundary.
  ! Discussion:
  ! owner( iProcFacesStart + i ) - Index of cell at this processor's boundary
  ! iOtherProc = iProcStart+i - Where the exchanged values of the boundary cells at 'other' processors are written.
  ! apr(i) - Matrix coefficients for cells that are on the processor boundary. Contributions from a cell that is in 'other' domain.
  do i=1,npro
    k = owner( iProcFacesStart + i )
    iOtherProc = iProcStart+i
    res( k ) = res( k ) - apr( i )*fi( iOtherProc )
  end do

  ! L1-norm of residual
  res0=sum(abs(res))
  call global_sum(res0)

  if(res0.lt.tol) then
    write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  PCG(IC0):  Solving for ',trim(chvarSolver(IFI)), &
    ', Initial residual = ',RES0,', Final residual = ',RES0,', No Iterations ',0
    return
  endif
!
! If ltest=true, print the norm 
!
  if(ltest) write(66,'(20x,a,1pe10.3)') 'res0 = ',res0

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
  call global_sum(sk)

!
! Calculate beta
!
  bet=sk/s0

!
! Calculate new search vector pk
!
  pk = zk + bet*pk

  call exchange( pk )
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

  ! Processor boundaries
  do i=1,npro
    k = owner( iProcFacesStart + i )
    iOtherProc = iProcStart+i
    zk( k ) = zk( k ) + apr( i )*pk( iOtherProc )
  end do

  ! Inner product
  pkapk=sum(pk*zk)
  call global_sum( pkapk )

  alf=sk/pkapk

  ! Update solution vector
  fi(1:numCells) = fi(1:numCells) + alf*pk

  ! Update residual vector
  res = res - alf*zk

  ! L1-norm of residual
  resl = sum(abs(res))
  call global_sum( resl )

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
  if ( myid.eq.0 ) write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') &
  'PCG(Jacobi):  Solving for ',trim(chvarSolver(IFI)), &
  ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L

end subroutine