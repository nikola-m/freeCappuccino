module types
  integer, parameter :: dp = kind(1.0d0)
end module


!***********************************************************************
!
subroutine bicgstab(a,ja,ioffset,diag,nnz,numCells,numTotal,su,fi)
!
!***********************************************************************
!
! BiCGStab for sparse matrices in CSR format.   
!
!***********************************************************************
!
  use types, only:dp
  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: nnz,numCells,numTotal,ja(nnz),ioffset(numCells+1),diag(numCells)
  real(dp), intent(in) :: a(nnz)
  real(dp), intent(in) :: su(numCells)
  real(dp), intent(inout) :: fi(numTotal)

!
!     Local variables
!
  integer :: i, k, ns, l
  real(dp), dimension(numCells) :: res,reso,pk,uk,zk,vk,d
  real(dp) :: rsm, resmax, res0, resl, resor
  real(dp) :: alf, beto, gam, bet, om, ukreso
  real(dp) :: small

! max no. of iterations
  resmax = 1e-13
  small  = 1e-30

!
! Calculate initial residual vector and the norm
!
  
  res = 0.0_dp

  do i=1,numCells
    res(i) = su(i) 
    do k = ioffset(i),ioffset(i+1)-1
      res(i) = res(i) -  a(k) * fi(ja(k)) 
    enddo
  enddo

  ! do i=1,noc
  !   res(ijl(i)) = res(ijl(i)) - ar(i)*fi(ijr(i))
  !   res(ijr(i)) = res(ijr(i)) - al(i)*fi(ijl(i))
  ! end do

  ! L^1-norm of residual
  res0=sum(abs(res))

!
! Print the norm 
!
  write(6,'(20x,a,1pe10.3)') 'res0 = ',res0

!
! Calculate elements of diagonal preconditioning matrix
!
  do i=1,numCells
    d(i) = a( diag(i) )
    do k = ioffset(i), diag(i)-1
      do l = diag( ja(k) ), ioffset( ja( k )+1 )-1
        ! kolona u kojoj je trenutni dijagonalni element je i
        ! kada je pronadje izlazi sa poslednjom vrednosti l indeksa:
        if ( ja( l ) == i ) exit 
      end do
      d(i) = d(i) - a( k ) * d( ja( k )) * a( l )
    end do
    d(i) = 1.0_dp / d(i)
  enddo 


!
! Initialize working arrays and constants
!
  reso  = res    
  pk = 0.0_dp
  uk = 0.0_dp
  zk = 0.0_dp
  vk = 0.0_dp

  alf = 1.0_dp
  beto = 1.0_dp
  gam = 1.0_dp


!
! Start iterations
!
  ns=100
  do l=1,ns

!
! Calculate beta and omega
!
  bet = sum(res*reso)
  om = bet*gam/(alf*beto+small)
  beto = bet


!
! Calculate pk
!
  pk = res + om*(pk -alf*uk)


!
! Solve (M ZK = PK) - forward substitution
!
  do i=1,numCells
    zk(i) = pk(i)
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


!
! Matvec 1: Uk = A*pk
!
  do i=1,numCells
    uk(i) = 0.0_dp
    do k = ioffset(i),ioffset(i+1)-1
      uk(i) = uk(i) +  a(k) * zk(ja(k)) 
    enddo
  enddo


!
! Calculate scalar product uk*reso, and gamma
!
  ukreso = sum(uk*reso)
  gam = bet/ukreso

!
! Update 'fi' and calculate 'w' (overwrite 'res; - it is res-update)
!
  fi(1:numCells) = fi(1:numCells) + gam*zk
  res            = res            - gam*uk   ! <- W

!
! Solve (M Y = W); Y overwrites zk; forward substitution
!
  do i=1,numCells
    zk(i) = res(i)
    do k = ioffset(i), diag(i)-1
      zk(i) = zk(i) -  a( k ) * zk( ja( k ) )
    end do
    zk(i) = zk(i)*d(i)
  enddo

  zk = zk/(d+small)

!
! Backward substitution
!
  do i=numCells,1,-1
    do k = diag(i)+1, ioffset(i+1)-1
      zk(i) = zk(i) - a( k ) * zk( ja( k ) )
    end do
    zk(i) = zk(i)*d(i)
  enddo


!
! Matvec 2: v = A*y (vk = A*zk); vk = csrMatVec(a,zk)
!
  do i=1,numCells
    vk(i) = 0.0_dp
    do k = ioffset(i),ioffset(i+1)-1
      vk(i) = vk(i) +  a(k) * zk(ja(k)) 
    enddo
  enddo  


!
! Calculate alpha (alf)
!
  alf = sum(vk*res) / (sum(vk*vk)+small)


  ! Update solution vector
  fi(1:numCells) = fi(1:numCells) + alf*zk

  ! Update residual vector
  res = res - alf*vk

  ! L^1-norm of residual
  resl = sum(abs(res))

!
! Check convergence
!
  if(l.eq.1) resor = res0
  rsm = resl/(resor+small)
  write(6,'(19x,a,i4,a,1pe10.3,a,1pe10.3)')  ' iter = ',l,' resl = ',resl,' rsm = ',rsm
  if(rsm.lt.resmax) exit

!
! End of iteration loop
!
  end do

end subroutine



!***********************************************************************
!
subroutine iccg(a,ja,ioffset,diag,nnz,numCells,numTotal,su,fi)
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
  use types, only:dp
  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: nnz,numCells,numTotal,ja(nnz),ioffset(numCells+1),diag(numCells)
  real(dp), intent(in) :: a(nnz)
  real(dp), intent(in) :: su(numCells)
  real(dp), intent(inout) :: fi(numTotal)

!
! Local variables
!
  integer :: i, k, ns, l
  real(dp), dimension(numCells) :: pk,zk,d,res
  real(dp) :: rsm, resmax, res0, resl, resor
  real(dp) :: s0, sk, alf, bet, pkapk, small

  small  = 1e-30

! max no. of iterations
  resmax = 1e-13
!
! Initalize working arrays
!
  pk(:) = 0.0_dp
  zk(:) = 0.0_dp
  d(:) = 0.0_dp
  res(:) = 0.0_dp
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

  ! do i=1,noc
  !   res(ijl(i)) = res(ijl(i)) - ar(i)*fi(ijr(i))
  !   res(ijr(i)) = res(ijr(i)) - al(i)*fi(ijl(i))
  ! end do

  ! L^1-norm of residual
  res0=sum(abs(res))

!
! If ltest=true, print the norm 
!
  write(6,'(20x,a,1pe10.3)') 'res0 = ',res0
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
  ns=100
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
  ! do i=1,noc
  !   zk(ijl(i)) = zk(ijl(i)) + ar(i)*pk(ijr(i))
  !   zk(ijr(i)) = zk(ijr(i)) + al(i)*pk(ijl(i))
  ! end do

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
  if(l.eq.1) resor = res0
  rsm = resl/(resor+small)
  write(6,'(19x,a,i4,a,1pe10.3,a,1pe10.3)') ' iter = ',l,' resl = ',resl,' rsm = ',rsm
  if(rsm.lt.resmax) exit
!
! End of iteration loop
!
  end do


end subroutine

!***********************************************************************
!
subroutine dpcg(a,ja,ioffset,diag,nnz,numCells,numTotal,su,fi)
!
!***********************************************************************
!
!    This routine incorporates the diagonally preconditioned 
!    Conjugate Gradient solver for symmetric matrices in CSR sparse matrix format
!
!    Writen by nikola mirkov, 2016. nmirkov@vinca.rs
!
!***********************************************************************
!
  use types, only:dp
  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: nnz,numCells,numTotal,ja(nnz),ioffset(numCells+1),diag(numCells)
  real(dp), intent(in) :: a(nnz)
  real(dp), intent(in) :: su(numCells)
  real(dp), intent(inout) :: fi(numTotal)

!
! Local variables
!
  integer :: i, k, ns, l
  real(dp), dimension(numCells) :: pk,zk,res
  real(dp) :: rsm, resmax, res0, resl, resor, small
  real(dp) :: s0, sk, alf, bet, pkapk

  small  = 1e-30

! max no. of iterations
  resmax = 1e-13

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

  ! do i=1,noc
  !   res(ijl(i)) = res(ijl(i)) - ar(i)*fi(ijr(i))
  !   res(ijr(i)) = res(ijr(i)) - al(i)*fi(ijl(i))
  ! end do

  ! L^1-norm of residual
  res0=sum(abs(res))

!
! Print the norm 
!
  write(6,'(20x,a,1pe10.3)') 'res0 = ',res0

  s0=1.e20
!
! Start iterations
!
  ns=100
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
  ! do i=1,noc
  !   zk(ijl(i)) = zk(ijl(i)) + ar(i)*pk(ijr(i))
  !   zk(ijr(i)) = zk(ijr(i)) + al(i)*pk(ijl(i))
  ! end do

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
  if(l.eq.1) resor = res0
  rsm = resl/(resor+small)
  write(6,'(19x,a,i4,a,1pe10.3,a,1pe10.3)') ' iter = ',l,' resl = ',resl,' rsm = ',rsm
  if(rsm.lt.resmax) exit
!
! End of iteration loop
!
  end do

end subroutine



program test_sparse
! The program computes the solution to the system of linear
! equations with a square matrix A and
! right-hand side B, where A is the coefficient matrix.
!
! Test 1 - Real nonsymmetric 5x5 system :
!
!   6.80  -6.05  -0.45   8.32  -9.67
!  -2.11  -3.30   2.58   2.71  -5.14
!   5.66   5.36  -2.70   4.35  -7.26
!   5.97  -4.44   0.27  -7.17   6.08
!   8.23   1.08   9.04   2.14  -6.87
!
! and B is the right-hand side:
!
!   4.02
!   6.19 
!  -8.22
!  -7.57
!  -3.03
!
! Solution
!  -0.80
!  -0.70
!   0.59
!   1.32 
!   0.57
!
!
! Test 2 - symmetric positive definite 5x5 system!
!
!    3.14   0.17  -0.90   1.65  -0.72
!    0.17   0.79   0.83  -0.65   0.28
!   -0.90   0.83   4.53  -3.70   1.60
!    1.65  -0.65  -3.70   5.32  -1.37
!   -0.72   0.28   1.60  -1.37   1.98
!
!  and B is the right-hand side:
!
!   -7.29
!    9.25
!    5.99
!   -1.94
!   -8.30 
!
! Solution
!  -6.02
!  15.62
!   3.02
!   3.25
!  -8.78
!
  use types, only: dp
  implicit none

  integer :: i,ioffset(6),diag(5)
  integer :: ja(25)
  real(dp) :: a(25), b(5), x(5), xa(5)

!
! Test 1 - Real nonsymmetric 5x5 system
!
  write(*,'(a)') ' '
  write(*,'(a)') ' Test 1 - Real nonsymmetric 5x5 system.'
  write(*,'(a)') ' '

  a = (/ &
  6.80, -6.05,  -0.45,   8.32,  -9.67,  &
  -2.11,  -3.30,   2.58,   2.71, -5.14, &
  5.66,  5.36,  -2.70,   4.35,  -7.26,  &
  5.97,  -4.44,   0.27,  -7.17,   6.08, &
  8.23,   1.08,   9.04,   2.14,  -6.87 /)

  b = (/ 4.02, 6.19,-8.22,-7.57,-3.03 /)

  ja = (/ &
      1,2,3,4,5,&
      1,2,3,4,5,&
      1,2,3,4,5,&
      1,2,3,4,5,&
      1,2,3,4,5 &
      /)

  ioffset = (/ 1, 6, 11, 16, 21, 26 /)

  diag = (/ 1, 7, 13, 19, 25 /)

  xa = (/ -0.80, -0.70, 0.59, 1.32, 0.57 /)

  write(*,'(a)') ' '
  write(*,'(19x,a)') ' The ILU(0) preconditioned BiCGStab solver:'

  call BiCGStab(a,ja,ioffset,diag,25,5,5,b,x)

  write(*,'(a)') ' '
  write(*,'(19x,a,5x,a)') 'Solution: ','Analytic solution:'
  do i=1,5
    write(*,'(19x,f5.2,10x,f5.2)') x(i),xa(i)
  enddo

!
! Test 2 - symmetric positive definite 5x5 system
!
  write(*,'(a)') ' '
  write(*,'(a)') ' Test 2 - symmetric positive definite 5x5 system.'
  write(*,'(a)') ' '

  a = (/ &
    3.14,   0.17,  -0.90,   1.65,  -0.72, &
    0.17,   0.79,   0.83,  -0.65,   0.28, &
   -0.90,   0.83,   4.53,  -3.70,   1.60, &
    1.65,  -0.65,  -3.70,   5.32,  -1.37, &
   -0.72,   0.28,   1.60,  -1.37,   1.98 /)

  b = (/ -7.29, 9.25, 5.99, -1.94, -8.30 /)

  xa = (/  -6.02,15.62,3.02,3.25,-8.78 /)

  write(*,'(a)') ' '
  write(*,'(19x,a)') ' The Incomplete Cholesky preconditioned CG solver:'

  call ICCG(a,ja,ioffset,diag,25,5,5,b,x)

  write(*,'(a)') ' '
  write(*,'(19x,a,5x,a)') 'Solution: ','Analytic solution:'
  do i=1,5
    write(*,'(19x,f5.2,10x,f5.2)') x(i),xa(i)
  enddo

  write(*,'(a)') ' '
  write(*,'(19x,a)') ' The Diagonally preconditioned CG solver:'

  call DPCG(a,ja,ioffset,diag,25,5,5,b,x)

  write(*,'(a)') ' '
  write(*,'(19x,a,5x,a)') 'Solution: ','Analytic solution:'
  do i=1,5
    write(*,'(19x,f5.2,10x,f5.2)') x(i),xa(i)
  enddo

end program

