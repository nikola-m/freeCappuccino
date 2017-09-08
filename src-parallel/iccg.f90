!***********************************************************************
!
subroutine iccg(fi,ifi)
!
!***********************************************************************
!
!    This routine incorporates the incomplete Cholesky preconditioned 
!    Conjugate Gradient solver for symmetric matrices in CSR sparse matrix format
!
!    Writen by nikola mirkov, 2016. nmirkov@vinca.rs
!
!***********************************************************************
!
  use types 
  use parameters
  use geometry, only: numCells,numTotal,npro,noc,ijl,ijr,iProcFacesStart,iProcStart
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
  real(dp), dimension(numCells) :: zk,d
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
  d = 0.0_dp
  res = 0.0_dp
!
! Calculate initial residual vector and the norm
!

! res(:) = su(:) - sum( a(ioffset(1:numCells):ioffset(2:numCells+1)-1) * fi(ja(ioffset(1:numCells):ioffset(2:numCells+1)-1)) )

  do i=1,numCells
    res(i) = su(i) 
    do k = ioffset(i),ioffset(i+1)-1
      res(i) = res(i) -  a(k) * fi( ja(k) ) 
    enddo
  enddo

  ! Residual contribution for cells at OC-cut boundary
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

  ! if( res0.lt.sor(ifi) ) return
!
! If ltest=true, print the norm 
!
  if(ltest) write(6,'(20x,a,1pe10.3)') 'res0 = ',res0
!
! Calculate elements of diagonal preconditioning matrix
!
  do i=1,numCells
    d(i) = a( diag(i) )
    do k = ioffset(i), diag(i)-1
      d(i) = d(i) - a( k )**2 * d( ja( k )) 
    end do
    d(i) =  1.0_dp / d(i)
  enddo

  s0=1.e20
!
! Start iterations
!
  ns=nsw(ifi)
  do l=1,ns
!
! Solve for zk(ijk) -- forward substitution
!
  do i=1,numCells
    zk(i) = res(i)
    do k = ioffset(i), diag(i)-1
      zk(i) = zk(i) -  a( k ) * zk( ja( k ))
    end do
    zk(i) = zk(i)*d(i)
  enddo

  zk = zk/(d+small)     
!
! Backward substitution
!
  do i=numCells,1,-1
    do k = diag(i)+1, ioffset(i+1)-1
      zk(i) = zk(i) - a( k ) * zk( ja( k ))
    end do
    zk(i) = zk(i)*d(i)
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

  ! OC-cut boundaries
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
  call global_sum(pkapk)

  alf=sk/pkapk

  ! Update solution vector
  fi(1:numCells) = fi(1:numCells) + alf*pk

  ! Update residual vector
  res = res - alf*zk


  ! L^1-norm of residual
  resl = sum(abs(res))
  call global_sum(resl)
  
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

  ! MPI exchange
  call exchange( fi )

! Write linear solver report:
  if ( myid.eq.0 ) write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  PCG(IC0):  Solving for ',trim(chvarSolver(IFI)), &
  ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L

end subroutine
