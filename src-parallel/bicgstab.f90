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
  use geometry, only: numCells,numTotal,npro,ijl,ijr,iProcFacesStart,iProcStart
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
  integer :: i, k, ns, l, iOtherProc
  real(dp), dimension(numCells) :: reso,pk,uk,vk,d
  real(dp), dimension(numCells+npro) :: zk
  real(dp) :: rsm, resmax, res0, resl, tol
  real(dp) :: alf, beto, gam, bet, om, ukreso, svkres, svkvk

! residual tolerance
  resmax = sor(ifi)
  tol = 1e-13

!
! Calculate initial residual vector and the norm
!

! res(:) = su(:) - sum( a(ioffset(1:numCells):ioffset(2:numCells+1)-1) * fi(ja(ioffset(1:numCells):ioffset(2:numCells+1)-1)) )
  
  res = 0.0_dp

  ! Residual contribution for ineer cells
  do i=1,numCells
    res(i) = su(i) 
    do k = ioffset(i),ioffset(i+1)-1
      res(i) = res(i) -  a(k) * fi(ja(k)) 
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
        ! kolona u kojoj je trenutni dijagonalni element je "i"
        ! kada je pronadje izlazi sa poslednjom vrednosti "l" indeksa:
        if ( ja( l ) == i ) exit 
      end do
      d(i) = d(i) - a( k ) * d( ja( k )) * a( l )
    end do
    d(i) = 1.0_dp / d(i)
  enddo 

!
! Initialize working arrays and parameters
!
  reso  = res    
  pk = 0.0_dp
  uk = 0.0_dp
  zk = 0.0_dp
  vk = 0.0_dp

  ! Parameters
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
  call global_sum(bet)

  om = bet*gam/(alf*beto+small)
  beto = bet

!
! Calculate pk
!
  pk = res + om*(pk - alf*uk)


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

  call exchange ( zk )

!
! Matvec 1: Uk = A*pk
!
  do i=1,numCells
    uk(i) = 0.0_dp
    do k = ioffset(i),ioffset(i+1)-1
      uk(i) = uk(i) +  a(k) * zk( ja(k) ) 
    enddo
  enddo


  ! Processor boundaries
  do i=1,npro
    k = owner( iProcFacesStart + i )
    iOtherProc = iProcStart+i
    uk( k ) = uk( k ) + apr( i ) * zk( iOtherProc )
  end do

!
! Calculate scalar product uk*reso, and gamma
!
  ukreso = sum(uk*reso)
  call global_sum(ukreso)

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

  call exchange( zk )

!
! Matvec 2: v = A*y (vk = A*zk); vk = csrMatVec(a,zk)
!
  do i=1,numCells
    vk(i) = 0.0_dp
    do k = ioffset(i),ioffset(i+1)-1
      vk(i) = vk(i) +  a(k) * zk( ja(k) ) 
    enddo
  enddo 

  ! Processor boundaries
  do i=1,npro
    k = owner( iProcFacesStart + i )
    iOtherProc = iProcStart+i
    vk( k ) = vk( k ) + apr( i ) * zk( iOtherProc )
  end do


!
! Calculate alpha (alf)
!
  svkres = sum(vk*res)
  call global_sum( svkres )

  svkvk = sum(vk*vk)
  call global_sum( svkvk )

  alf = svkres / (svkvk+small)

  ! Update solution vector
  fi(1:numCells) = fi(1:numCells) + alf*zk

  ! Update residual vector
  res = res - alf*vk

  ! L1-norm of residual
  resl = sum(abs(res))
  call global_sum(resl)

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
  if ( myid.eq.0 ) write(6,'(3a,1PE10.3,a,1PE10.3,a,I0)') '  BiCGStab(ILU(0)):  Solving for ',trim(chvarSolver(IFI)), &
  ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L

end subroutine

