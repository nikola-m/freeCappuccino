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
  use indexes
  use sparse_matrix
  use title_mod

  implicit none
!
!***********************************************************************
!
  integer, intent(in) :: ifi
  real(prec), dimension(numTotal), intent(inout) :: fi 

!
! Local variables
!
  integer :: i, k, ns, l
  real(prec), dimension(numCells) :: pk,zk,d
  real(prec) :: rsm, resmax, res0, resl
  real(prec) :: s0, sk, alf, bet, pkapk

! max no. of iterations
  resmax = sor(ifi)
!
! Initalize working arrays
!
  pk(:) = 0.0d0
  zk(:) = 0.0d0
  d(:) = 0.0d0
  res(:) = 0.0d0
!
! Calculate initial residual vector and the norm
!
  ! do k=2,nkm
  !   do i=2,nim
  !     do j=2,njm
  !       ijk=lk(k)+li(i)+j

  !       res(ijk)=su(ijk)-(ae(ijk)*fi(ijk+nj) +aw(ijk)*fi(ijk-nj) &
  !                        +an(ijk)*fi(ijk+1)  +as(ijk)*fi(ijk-1) &
  !                        +at(ijk)*fi(ijk+nij)+ab(ijk)*fi(ijk-nij)) &
  !                       -ap(ijk)*fi(ijk)

  !     end do
  !   end do
  ! end do

res(1:numCells)=su(1:numCells)-sum(a(ioffset(1:numCells):ioffset(2:numCells+1)-1)*fi(ja(ioffset(1:numCells):ioffset(2:numCells+1)-1))) !- a( diag(1:numCells) )*fi(1:numCells)

!..or

do i=1,numCells
  res(i) = su(i) - sum(a( ioffset(i):ioffset(i+1)-1 )*fi(ja(ioffset(i):ioffset(i+1)-1))) !- a( diag(i) )*fi(i)
enddo

!..or

do i=1,numCells
  res(i) = su(i) 
  do k = ioffset(i),ioffset(i+1)-1
    res(i) = res(i) - a( k )*fi(ja(k)) 
  enddo
enddo

!.. or C variant
! for(int i=0; i < numCells, i++){
!   res[i] = su[i]
!   for(int k=ioffset[i]; i<ioffset[i]; i++){
!     res[i] = res[i] - (a[k]*fi[ja[k]])
!   }
! }



! ..or
res(owner(:)) = res(owner(:)) + a( icell_jcell_csr_value_index(:) )*fi(neighbour(:))
res(neighbour(:)) = res(neighbour(:)) + a( jcell_icell_csr_value_index(:) )*fi(owner(:))

res(:) = su(:) - res(:)- a(diag(:))*fi(1:numCells)

res(ijl(:)) = res(ijl(:)) - ar(:)*fi(ijr(:))
res(ijr(:)) = res(ijr(:)) - al(:)*fi(ijl(:))

!..or

  ! ! res_i = su_i - sum_j(a_{i,j}*phi_j)
  ! do i = 1,numInnerFaces
  !   ijp = owner(i)
  !   ijn = neighbour(i)

  !     k = icell_jcell_csr_value_index(i)
  !     res(ijp) = res(ijp) + a(k)*fi(ijn)

  !     k = jcell_icell_csr_value_index(i)
  !     res(ijn) = res(ijn) + a(k)*fi(ijp)

  ! enddo
  ! do ijp=1,numCells
  !     res(ijp) = su(ijp) - res(ijp)- a(diag(ijp))*fi(ijp)
  ! enddo

  ! do i=1,noc
  !   res(ijl(i)) = res(ijl(i)) - ar(i)*fi(ijr(i))
  !   res(ijr(i)) = res(ijr(i)) - al(i)*fi(ijl(i))
  ! end do





  ! L^1-norm of residual
  res0=sum(abs(res(:)))

!
! If ltest=true, print the norm 
!
  if(ltest) write(66,'(20x,a,1pe10.3)') 'res0 = ',res0
!
!.....calculate elements of diagonal preconditioning matrix
!
      ! do k=2,nkm
      !   do i=2,nim
      !     do j=2,njm
      !       ijk=lk(k)+li(i)+j
      !       d(ijk)=1./(ap(ijk) &
      !         -aw(ijk)**2*d(ijk-nj) &
      !         -as(ijk)**2*d(ijk-1)  &
      !         -ab(ijk)**2*d(ijk-nij))
      !     end do
      !   end do
      ! end do

  !matvec - samo sa elementima koji su LEVO od dijagonale: 
  do i=1,numCells
    d(i) = a( diag(i) )
    do k = ioffset(i), diag(i)-1
      d(i) = d(i) - a( k )**2 * d( ja( k )) 
    end do
    d(i) =  1. / d(i)
  enddo

!..or

  ! !matvec - samo sa elementima koji su LEVO od dijagonale: 
  ! do i=1,numCells
  !   d(i) = 0.0d0
  !   row_elements_loop: do k = ioffset(i), ioffset(i+1)-1
  !     if ( k < diag(i)) then !ili ja(k) < icell, ali ovako nije potrebno jer nam je 'ja' u rastucem nizu
  !       d(i) = d(i) - a( k )**2 * d( ja( k )) 
  !     elseif ( k = diag(i)) then
  !       d(i) = d(i) + 1./a( diag(i))
  !     else
  !       exit row_elements_loop ! ili cycle ako bi gore koristili ja(k), ali ovako nije potrebno jer nam je 'ja' u rastucem nizu
  !     endif 
  !   end do
  ! enddo


  s0=1.e20
!
! Start iterations
!
  ns=nsw(ifi)
  do l=1,ns
!
!.....solve for zk(ijk) -- forward substitution
!
      ! do k=2,nkm
      !   do i=2,nim
      !     do j=2,njm
      !       ijk=lk(k)+li(i)+j
      !       zk(ijk)=(res(ijk)-aw(ijk)*zk(ijk-nj)-as(ijk)*zk(ijk-1)-  &
      !               ab(ijk)*zk(ijk-nij))*d(ijk)
      !     end do
      !   end do
      ! end do

  !matvec - samo sa elementima koji su LEVO od dijagonale: 
  do i=1,numCells
    zk(i) = res(i)
    do k = ioffset(i), diag(i)-1
      zk(i) = zk(i) -  a( k )*zk( ja( k ))
    end do
    zk(i) = zk(i)*d(i)
  enddo

!..or

  ! !matvec - samo sa elementima koji su LEVO od dijagonale: 
  ! do i=1,numCells
  !   vec(i) = res(i)
  !   row_elements_loop: do k = ioffset(i), ioffset(i+1)-1
  !     if ( k < diag(i)) then
  !       vec(i) = vec(i) +  a( k )*zk( ja( k ))
  !     else
  !       exit row_elements_loop
  !     endif 
  !   end do
  !   zk(i) = zk(i)*d(i)
  ! enddo


      ! do k=2,nkm
      !   do i=2,nim
      !     do j=2,njm
      !       ijk=lk(k)+li(i)+j
      !       zk(ijk)=zk(ijk)/(d(ijk)+small)
      !     end do
      !   end do
      ! end do

  zk(:) = zk(:)/(d(:)+small)     
!
!..... backward substitution; calculate scalar product sk
!
      ! sk=0.0d0
      ! do k=nkm,2,-1
      !   do i=nim,2,-1
      !     do j=njm,2,-1
      !       ijk=lk(k)+li(i)+j
      !       zk(ijk)=(zk(ijk)-(ae(ijk)*zk(ijk+nj)+an(ijk)*zk(ijk+1)+ &
      !                at(ijk)*zk(ijk+nij)))*d(ijk)
      !       sk=sk+res(ijk)*zk(ijk)
      !     end do
      !   end do
      ! end do

  !matvec - samo sa elementima koji su DESNO od dijagonale: 
  do i=numCells,1,-1
    zk(i) = zk(i)
    do k = diag(i)+1, ioffset(i+1)-1
      zk(i) = zk(i) - ( a( k ) * zk( ja( k )) )
    end do
    zk(i) = zk(i)*d(i)
  enddo

!..or

  ! !matvec - samo sa elementima koji su DESNO od dijagonale: 
  ! do i=numCells,1,-1
  !   zk(i) = zk(i)
  !   row_elements_loop: do k = ioffset(i), ioffset(i+1)-1
  !     if ( k .le. diag(i)) then !ili ja(k) < icell, ali ovako nije potrebno jer nam je ja u rastucem nizu
  !       cycle row_elements_loop
  !     else 
  !       zk(i) = zk(i) - ( a( k ) * zk( ja( k )) )
  !     endif 
  !   end do
  !   zk(i) = zk(i)*d(i)
  ! enddo

  sk = sum(res(:)*zk(:))
!
! Calculate beta
!
  bet=sk/s0
!
!.....calculate new search vector pk
!
      ! do k=2,nkm
      !   do i=2,nim
      !     do j=2,njm
      !       ijk=lk(k)+li(i)+j
      !       pk(ijk)=zk(ijk)+bet*pk(ijk)
      !     end do
      !   end do
      ! end do

  pk(:) = zk(:) + bet*pk(:)
!
!.... calculate scalar product (pk.a pk) and alpha (overwrite zk)
!
      ! pkapk=0.0d0
      ! do k=2,nkm
      !   do i=2,nim
      !     do j=2,njm
      !       ijk=lk(k)+li(i)+j
      !       zk(ijk)=ap(ijk)*pk(ijk)+ae(ijk)*pk(ijk+nj)+                &
      !         aw(ijk)*pk(ijk-nj)+an(ijk)*pk(ijk+1)+as(ijk)*pk(ijk-1)+  &
      !         at(ijk)*pk(ijk+nij)+ab(ijk)*pk(ijk-nij)
      !       pkapk=pkapk+pk(ijk)*zk(ijk)
      !     end do
      !   end do
      ! end do

! 1)
do i=1,numCells
  zk(i) = 0.0d0 
  do k = ioffset(i),ioffset(i+1)-1
    zk(i) = zk(i) + a( k )*pk(ja(k)) 
  enddo
enddo

!  !..or 2)
!   do i=1,numCells
!    zk(i) = sum( a( ioffset(i) : ioffset(i+1)-1 ) * pk( ja( ioffset(i) : ioffset(i+1)-1 )) ) 
!   enddo

! !..or 3)
!   zk(1:numCells) = sum( a(ioffset(1:numCells):ioffset(2:numCells+1)-1)*pk(ja(ioffset(1:numCells):ioffset(2:numCells+1)-1)) )
 

  ! block-boundaries
  do i=1,noc
    zk(ijl(i)) = zk(ijl(i)) + ar(i)*pk(ijr(i))
    zk(ijr(i)) = zk(ijr(i)) + al(i)*pk(ijl(i))
  end do

  ! zk(ijl(:)) = zk(ijl(:)) + ar(i)*pk(ijr(:))
  ! zk(ijr(:)) = zk(ijr(:)) + al(i)*pk(ijl(:))






! ! ..or 4) in face form:
!   do i = 1,numInnerFaces
!     ijp = owner(i)
!     ijn = neighbour(i)
!       k = icell_jcell_csr_value_index(i)
!       zk(ijp) = zk(ijp) + a(k)*pk(ijn)

!       k = jcell_icell_csr_value_index(i)
!       zk(ijn) = zk(ijn) + a(k)*pk(ijp)
!   enddo

!   do ijp=1,numCells
!       zk(ijp) = zk(ijp) + a(diag(ijp))*pk(ijp)
!   enddo

!   do i=1,noc
!     zk(ijl(i)) = zk(ijl(i)) + ar(i)*pk(ijr(i))
!     zk(ijr(i)) = zk(ijr(i)) + al(i)*pk(ijl(i))
!   end do


! !..or 5) in face form vectorised
!   zk(owner(:)) = zk(owner(:))         + ( a(icell_jcell_csr_value_index(:))*pk(neighbour(:)) )
!   zk(neighbour(:)) = zk(neighbour(:)) + ( a(jcell_icell_csr_value_index(:))*pk(owner(:))     )
!   zk(:) = zk(:) + a(diag(:))*pk(:)
!   ! block-boundaries
!   zk(ijl(:)) = zk(ijl(:)) + ar(i)*pk(ijr(:))
!   zk(ijr(:)) = zk(ijr(:)) + al(i)*pk(ijl(:))





  ! Inner product
  pkapk=sum(pk(:)*zk(:))

  alf=sk/pkapk

  ! Update solution vector
  fi(:) = fi(:) + alf*pk(:)

  ! Update residual vector
  res(:) = res(:) - alf*zk(:)

  ! L^1-norm of residual
  resl = sum(abs(res(:)))

  s0=sk
!
! Check convergence
!
  if(l.eq.1) resor(ifi) = res0
  rsm = resl/(resor(ifi)+small)
  if(ltest) write(66,'(19x,3a,i4,a,1pe10.3,a,1pe10.3)') ' fi=',chvar(ifi),' sweep = ',l,' resl = ',resl,' rsm = ',rsm
  if(rsm.lt.resmax) exit
!
! End of iteration loop
!
  end do

! Write linear solver report:
  write(66,'(3a,1PE10.3,a,1PE10.3,a,I0)') &
  'PCG(IC0):  Solving for ',trim(chvarSolver(IFI)), &
  ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L

end subroutine
