subroutine bicgstab(fi,ifi)
!
!***********************************************************************
!
! BiCGStab for sparse matrices in CSR format.   
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
  real(dp), dimension(numTotal) :: fi 

!
!     Local variables
!
  integer :: i, k, ns, l
  real(dp), dimension(numCells) :: reso,pk,uk,zk,vk,d
  real(dp) :: rsm, resmax, res0, resl
  real(dp) :: alf, beto, gam, bet, om, ukreso

! residual tolerance
  resmax = sor(ifi)

!
! Calculate initial residual vector and the norm
!

! res(:) = su(:) - sum( a(ioffset(1:numCells):ioffset(2:numCells+1)-1) * fi(ja(ioffset(1:numCells):ioffset(2:numCells+1)-1)) )
  
  res(:) = 0.0_dp

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
  res0=sum(abs(res(:)))
    if(res0.lt.small) then
      write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  BiCGStab(ILU(0)):  Solving for ',trim(chvarSolver(IFI)), &
      ', Initial residual = ',RES0,', Final residual = ',RES0,', No Iterations ',0
      return
    endif
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
  ns=nsw(ifi)
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
  ukreso = sum(uk(:)*reso(:))
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
  if(l.eq.1) resor(ifi) = res0
  rsm = resl/(resor(ifi)+small)
  if(ltest) write(6,'(19x,3a,i4,a,1pe10.3,a,1pe10.3)') ' fi=',chvar(ifi),' sweep = ',l,' resl = ',resl,' rsm = ',rsm
  if(rsm.lt.resmax) exit

!
! End of iteration loop
!
  end do

! Write linear solver report:
  write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  BiCGStab(ILU(0)):  Solving for ',trim(chvarSolver(IFI)), &
  ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L

end subroutine

