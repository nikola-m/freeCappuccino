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
  use geometry, only: numCells,numTotal,ijl,ijr
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
  real(dp) :: rsm, resmax, res0, resl


! max no. of iterations
  resmax = sor(ifi)

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
    fi(i) = fi(i) + res(i)/a(diag(i))   
  enddo

  do i=1,noc
    res(ijl(i)) = res(ijl(i)) - ar(i)*fi(ijr(i))
    res(ijr(i)) = res(ijr(i)) - al(i)*fi(ijl(i))
  end do

! L^1-norm of residual
  if(l.eq.1)  res0=sum(abs(res))

!
! If ltest=true, print the norm 
!
  if(ltest.and.l.eq.1) write(6,'(20x,a,1pe10.3)') 'res0 = ',res0 

!
! Update variable
!
  ! do i=1,numcells
  !   fi(i) = fi(i) + res(i)/a(diag(i))
  ! enddo

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

 

  ! L^1-norm of residual
  resl = sum(abs(res))

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
  write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  Gauss-Seidel:  Solving for ',trim(chvarSolver(IFI)), &
  ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L  

end subroutine