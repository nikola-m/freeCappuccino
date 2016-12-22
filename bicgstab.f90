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

! max no. of iterations
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

 
! !
! !.....CALCULATE ELEMENTS OF PRECONDITIONING MATRIX DIAGONAL
! !
!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             D(IJK)=1./(AP(IJK) - AW(IJK)*D(IJK-NJ)*AE(IJK-NJ) &
!                    - AS(IJK)*D(IJK-1)*AN(IJK-1)               &
!                    - AB(IJK)*D(IJK-NIJ)*AT(IJK-NIJ)) 
!           END DO
!         END DO
!       END DO



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

! !
! !..... CALCULATE BETA AND OMEGA
! !
!       BET=0.0_dp
!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             BET=BET+RES(IJK)*RESO(IJK)
!           END DO
!         END DO
!       END DO
!       OM=BET*GAM/(ALF*BETO+SMALL)
!       BETO=BET

!
! Calculate beta and omega
!
  bet = sum(res*reso)
  om = bet*gam/(alf*beto+small)
  beto = bet


! !
! !..... CALCULATE PK
! !
!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             PK(IJK)=RES(IJK)+OM*(PK(IJK)-ALF*UK(IJK))
!           END DO
!         END DO
!       END DO

!
! Calculate pk
!
  pk = res + om*(pk -alf*uk)

! !
! !.....SOLVE (M ZK = PK) - FORWARD SUBSTITUTION
! !
!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             ZK(IJK)=(PK(IJK)-AW(IJK)*ZK(IJK-NJ)  &
!                             -AS(IJK)*ZK(IJK-1)   &
!                             -AB(IJK)*ZK(IJK-NIJ))*D(IJK)
!           END DO
!         END DO
!       END DO

!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             ZK(IJK)=ZK(IJK)/(D(IJK)+SMALL)
!           END DO
!         END DO
!       END DO

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

! !
! !..... BACKWARD SUBSTITUTION
! !
!       DO K=NKM,2,-1
!         DO I=NIM,2,-1
!           DO J=NJM,2,-1
!             IJK=LK(K)+LI(I)+J
!             ZK(IJK)=(ZK(IJK)-AE(IJK)*ZK(IJK+NJ)  &
!                             -AN(IJK)*ZK(IJK+1)   &
!                             -AT(IJK)*ZK(IJK+NIJ))*D(IJK)
!           END DO
!         END DO
!       END DO

!
! Backward substitution
!
  do i=numCells,1,-1
    do k = diag(i)+1, ioffset(i+1)-1
      zk(i) = zk(i) - a( k ) * zk( ja( k ))
    end do
    zk(i) = zk(i)*d(i)
  enddo


! !
! !.....CALCULATE UK = A.PK
! !
!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             UK(IJK)=AP(IJK)*ZK(IJK)+AE(IJK)*ZK(IJK+NJ) +AW(IJK)*ZK(IJK-NJ) &
!                                    +AN(IJK)*ZK(IJK+1)  +AS(IJK)*ZK(IJK-1)  &
!                                    +AT(IJK)*ZK(IJK+NIJ)+AB(IJK)*ZK(IJK-NIJ)
!           END DO
!         END DO
!       END DO

!
! Matvec 1: Uk = A*pk
!
  do i=1,numCells
    uk(i) = 0.0_dp
    do k = ioffset(i),ioffset(i+1)-1
      uk(i) = uk(i) +  a(k) * zk(ja(k)) 
    enddo
  enddo

! !
! !..... CALCULATE SCALAR PRODUCT UK.RESO AND GAMMA
! !
!       UKRESO=0.0_dp
!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             UKRESO=UKRESO+UK(IJK)*RESO(IJK)
!           END DO
!         END DO
!       END DO
!       GAM=BET/UKRESO

!
! Calculate scalar product uk*reso, and gamma
!
  ukreso = sum(uk(:)*reso(:))
  gam = bet/ukreso

! !
! !.....UPDATE (FI) AND CALCULATE W (OVERWRITE RES - IT IS RES-UPDATE)
! !
!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             FI(IJK)=FI(IJK)+GAM*ZK(IJK)   
!             RES(IJK)=RES(IJK)-GAM*UK(IJK) !W
!           END DO
!         END DO
!       END DO

!
! Update 'fi' and calculate 'w' (overwrite 'res; - it is res-update)
!
  fi(1:numCells) = fi(1:numCells) + gam*zk
  res            = res            - gam*uk   ! <- W

! !
! !.....SOLVE (M Y = W); Y OVERWRITES ZK; FORWARD SUBSTITUTION
! !
!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             ZK(IJK)=(RES(IJK)-AW(IJK)*ZK(IJK-NJ) &
!                              -AS(IJK)*ZK(IJK-1)  &
!                              -AB(IJK)*ZK(IJK-NIJ))*D(IJK)
!            END DO
!          END DO
!       END DO

!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             ZK(IJK)=ZK(IJK)/(D(IJK)+SMALL)
!           END DO
!         END DO
!       END DO

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

! !
! !.....BACKWARD SUBSTITUTION
! !
!       DO K=NKM,2,-1
!         DO I=NIM,2,-1
!           DO J=NJM,2,-1
!             IJK=LK(K)+LI(I)+J
!             ZK(IJK)=(ZK(IJK)-AE(IJK)*ZK(IJK+NJ)  &
!                             -AN(IJK)*ZK(IJK+1)   &
!                             -AT(IJK)*ZK(IJK+NIJ))*D(IJK)
!           END DO
!         END DO
!       END DO

!
! Backward substitution
!
  do i=numCells,1,-1
    do k = diag(i)+1, ioffset(i+1)-1
      zk(i) = zk(i) - a( k ) * zk( ja( k ) )
    end do
    zk(i) = zk(i)*d(i)
  enddo

! !
! !.....CALCULATE V = A.Y (VK = A.ZK)
! !
!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             VK(IJK)=AP(IJK)*ZK(IJK)   +AE(IJK)*ZK(IJK+NJ)+  &
!                     AW(IJK)*ZK(IJK-NJ)+AN(IJK)*ZK(IJK+1)+   &
!                     AS(IJK)*ZK(IJK-1) +AT(IJK)*ZK(IJK+NIJ)+ &
!                     AB(IJK)*ZK(IJK-NIJ)
!           END DO
!         END DO
!       END DO

!
! Matvec 2: v = A*y (vk = A*zk); vk = csrMatVec(a,zk)
!
  do i=1,numCells
    vk(i) = 0.0_dp
    do k = ioffset(i),ioffset(i+1)-1
      vk(i) = vk(i) +  a(k) * zk(ja(k)) 
    enddo
  enddo  

! !
! !..... CALCULATE ALPHA (ALF)
! !
!       VRES=0.0_dp
!       VV=0.0_dp
!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             VRES=VRES+VK(IJK)*RES(IJK)
!             VV=VV+VK(IJK)*VK(IJK)
!           END DO
!         END DO
!       END DO

!       ALF=VRES/(VV+SMALL)
!
! Calculate alpha (alf)
!
  alf = sum(vk*res) / (sum(vk*vk)+small)

! !
! !.....UPDATE VARIABLE (FI) AND RESIDUAL (RES) VECTORS
! !
!       RESL=0.0_dp
!       DO K=2,NKM
!         DO I=2,NIM
!           DO J=2,NJM
!             IJK=LK(K)+LI(I)+J
!             FI(IJK)=FI(IJK)+ALF*ZK(IJK)
!             RES(IJK)=RES(IJK)-ALF*VK(IJK)
!             RESL=RESL+ABS(RES(IJK))
!           END DO
!         END DO
!       END DO

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

