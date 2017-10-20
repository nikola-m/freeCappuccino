!***********************************************************************
!
subroutine get_rAU_x_UEqnH()
!
!***********************************************************************
!
!     Assemble A(U), and H(U) excluding pressure gradient
!     Update U according to U = 1/A(U) * H(U)
!     Peric denotes such U by Utilde
!
!***********************************************************************
!
  use types
  use parameters
  use geometry
  use variables
  use sparse_matrix
  use gradients
  use hcoef, only: h

  implicit none
!
!***********************************************************************
!

! Local variables
  integer :: i, k, ijp, ijn, inp
  real(dp) :: apotime,sut,svt,swt
  real(dp) :: heat
  real(dp) :: off_diagonal_terms

  su = 0.0_dp
  sv = 0.0_dp
  sw = 0.0_dp

  ! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
  do inp=1,numCells

        !=====================================
        !.....BUOYANCY SOURCE TERMS
        !=====================================
        if(lcal(ien).and.lbuoy) then
        !----------------------------------------------
        !........[Boussinesq-ova aproximacija: ]
        !----------------------------------------------
          heat=0.0_dp
          if(boussinesq) then
           heat=beta*densit*(t(inp)-tref)*vol(inp)
          else
            heat=(densit-den(inp))*vol(inp)
          endif
        !----------------------------------------------
          su(inp)=su(inp)-gravx*heat
          sv(inp)=sv(inp)-gravy*heat
          sw(inp)=sw(inp)-gravz*heat
        endif


        !=====================================
        !.....UNSTEADY TERM
        !=====================================
        if(bdf) then
        !----------------------------------------------
        !    Three Level Implicit Time Integration Method:
        !    in case that BTIME=0. --> Implicit Euler
        !----------------------------------------------
              apotime=den(inp)*vol(inp)/timestep

              sut = apotime*((1+btime)*uo(inp))
              svt = apotime*((1+btime)*vo(inp))
              swt = apotime*((1+btime)*wo(inp))
            
              if (btime > 0.99) then ! bdf2 scheme btime=1.
                sut = sut - apotime*(0.5*btime*uoo(inp))
                svt = svt - apotime*(0.5*btime*voo(inp))
                swt = swt - apotime*(0.5*btime*woo(inp))
              endif

              su(inp) = su(inp) + sut
              sv(inp) = sv(inp) + svt
              sw(inp) = sw(inp) + swt

              ! spu(inp) = spu(inp) + apotime*(1+0.5*btime)
              ! spv(inp) = spv(inp) + apotime*(1+0.5*btime)
              ! sp(inp)  = sp(inp)  + apotime*(1+0.5*btime)

        endif

  end do

  !
  ! U component
  !
  if(cn) then

      do i = 1,numInnerFaces
          ijp = owner(i)
          ijn = neighbour(i)

          k = icell_jcell_csr_value_index(i)
          su(ijp) = su(ijp) - h(k)*uo(ijn)

          k = jcell_icell_csr_value_index(i)
          su(ijn) = su(ijn) - h(k)*uo(ijp)

      enddo

      do ijp=1,numCells
          apotime=den(ijp)*vol(ijp)/timestep
          off_diagonal_terms = sum( h( ioffset(ijp) :  ioffset(ijp+1)-1 ) ) - h(diag(ijp))
          su(ijp) = su(ijp) + (apotime + off_diagonal_terms)*uo(ijp)
      enddo

  endif

  ! Assemble H(U) = sum_j {a_j*U_pj}, j - runs trough neighbour indices
  do i = 1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      k = icell_jcell_csr_value_index(i)
      su(ijp) = su(ijp) + h(k)*u(ijn)

      k = jcell_icell_csr_value_index(i)
      su(ijn) = su(ijn) + h(k)*u(ijp)

  enddo

  !
  ! V component
  !
  if(cn) then

      do i = 1,numInnerFaces
          ijp = owner(i)
          ijn = neighbour(i)

          k = icell_jcell_csr_value_index(i)
          sv(ijp) = sv(ijp) - h(k)*vo(ijn)

          k = jcell_icell_csr_value_index(i)
          sv(ijn) = sv(ijn) - h(k)*vo(ijp)

      enddo

      do ijp=1,numCells
          apotime=den(ijp)*vol(ijp)/timestep
          off_diagonal_terms = sum( h( ioffset(ijp) :  ioffset(ijp+1)-1 ) ) - h(diag(ijp))
          sv(ijp) = sv(ijp) + (apotime + off_diagonal_terms)*vo(ijp)
      enddo

  endif

  ! Assemble H(V) = sum_j {a_j*V_pj}, j - runs trough neighbour indices
  do i = 1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      k = icell_jcell_csr_value_index(i)
      sv(ijp) = sv(ijp) + h(k)*v(ijn)

      k = jcell_icell_csr_value_index(i)
      sv(ijn) = sv(ijn) + h(k)*v(ijp)

  enddo

  !
  ! W component
  !
  if(cn) then

      do i = 1,numInnerFaces
          ijp = owner(i)
          ijn = neighbour(i)

          k = icell_jcell_csr_value_index(i)
          sw(ijp) = sw(ijp) - h(k)*wo(ijn)

          k = jcell_icell_csr_value_index(i)
          sw(ijn) = sw(ijn) - h(k)*wo(ijp)

      enddo

      do ijp=1,numCells
          apotime=den(ijp)*vol(ijp)/timestep
          off_diagonal_terms = sum( h( ioffset(ijp) :  ioffset(ijp+1)-1 ) ) - h(diag(ijp))
          sw(ijp) = sw(ijp) + (apotime + off_diagonal_terms)*wo(ijp)
      enddo

  endif

  ! Assemble H(W) = sum_j {a_j*W_pj}, j - runs trough neighbour indices
  do i = 1,numInnerFaces
      ijp = owner(i)
      ijn = neighbour(i)

      k = icell_jcell_csr_value_index(i)
      sw(ijp) = sw(ijp) + h(k)*w(ijn)

      k = jcell_icell_csr_value_index(i)
      sw(ijn) = sw(ijn) + h(k)*w(ijp)

  enddo



  ! Finally U = rAU*UEqnH()
  u(1:numCells) = apu*su
  v(1:numCells) = apv*sv
  w(1:numCells) = apw*sw

end subroutine
