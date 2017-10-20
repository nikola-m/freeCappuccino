!***********************************************************************
!
subroutine GaussSeidel(fi,ifi)
!
!***********************************************************************
!
!    This routine incorporates the  
!    Gauss-Seidel solver for sparse matrices in CSR sparse matrix format
!
!    Writen by nikola mirkov, 2016. nmirkov@vinca.rs
!
!***********************************************************************
!
  use types 
  use parameters
  use geometry, only: numCells,numTotal!,ijl,ijr
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
  integer :: i, k, ns, l, itr_used
  real(dp) :: rsm, resmax, res0, resl, tol


! residual tolerance
  resmax = sor(ifi)
  tol = 1e-13

  itr_used = 0

!
! Start iterations
!
  ns=nsw(ifi)
  do l=1,ns

!
! Calculate initial residual vector and the norm
!
  do i=1,numCells
    res(i) = su(i) 
    do k = ioffset(i),ioffset(i+1)-1
      res(i) = res(i) -  a(k) * fi(ja(k)) 
    enddo
    fi(i) = fi(i) + res(i)/(a(diag(i))+small)   
  enddo

  ! do i=1,noc
  !   res(ijl(i)) = res(ijl(i)) - ar(i)*fi(ijr(i))
  !   res(ijr(i)) = res(ijr(i)) - al(i)*fi(ijl(i))
  ! end do

! L1-norm of residual
  if(l.eq.1)  then
    res0=sum(abs(res))

    if(res0.lt.tol) then
      write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  Gauss-Seidel:  Solving for ',trim(chvarSolver(ifi)), &
      ', Initial residual = ',res0,', Final residual = ',res0,', No Iterations ',0
      return
    endif  

  endif

  if(ltest.and.l.eq.1) write(6,'(20x,a,1pe10.3)') 'res0 = ',res0  

  ! L1-norm of residual
  resl = sum(abs(res))

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
  write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  Gauss-Seidel:  Solving for ',trim(chvarSolver(ifi)), &
  ', Initial residual = ',res0,', Final residual = ',resl,', No Iterations ',itr_used 

end subroutine