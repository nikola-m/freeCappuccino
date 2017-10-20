!***********************************************************************
!
subroutine Jacobi(fi,ifi)
!
!***********************************************************************
!
!    This routine incorporates the  
!    Jacobi solver for sparse matrices in CSR sparse matrix format
!
!    Writen by nikola mirkov, 2016. nmirkov@vinca.rs
!
!***********************************************************************
!
  use types 
  use parameters
  use geometry, only: numCells,numTotal,ijl,ijr,noc,npro,iProcFacesStart,iProcStart
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
  integer :: i, k, ns, l, itr_used, iOtherProc
  real(dp) :: rsm, resmax, res0, resl, tol


! residual tolerance
  tol = 1e-13
  resmax = sor(ifi)

  itr_used = 0

!
! Calculate initial residual vector and the norm
!
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

  do i=1,npro
    k = owner( iProcFacesStart + i )
    iOtherProc = iProcStart+i
    res( k ) = res( k ) - apr( i )*fi( iOtherProc )
  end do


  ! L1-norm of residual
  res0=sum(abs(res))

  call global_sum(res0)

  if(res0.lt.tol) then
    if (myid .eq. 0) write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  Jacobi:  Solving for ',trim(chvarSolver(ifi)), &
    ', Initial residual = ',res0,', Final residual = ',res0,', No Iterations ',0
    return
  endif

  if(ltest) write(6,'(20x,a,1pe10.3)') 'res0 = ',res0 

!
! Start iterations
!
  ns=nsw(ifi)

  do l=1,ns

!
! Update solution vector
!
  do i=1,numCells
    fi(i) = fi(i) + res(i)/(a(diag(i))+small)   
  enddo

  call exchange(fi)

!
! Update residual vector
!
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

  do i=1,npro
    k = owner( iProcFacesStart + i )

    iOtherProc = iProcStart+i
    res( k ) = res( k ) - apr( i )*fi( iOtherProc )
  end do

  ! L1-norm of residual
  resl = sum(abs(res))

  call global_sum( resl )

  itr_used = itr_used + 1
  
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
  '  Jacobi:  Solving for ',trim(chvarSolver(ifi)), &
  ', Initial residual = ',res0,', Final residual = ',resl,', No Iterations ',itr_used 

end subroutine