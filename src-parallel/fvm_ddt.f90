subroutine fvm_ddt_vector_field
!  
!******************************************************************************
!
!     Finds source and matrix coefficients representing FVM discretization 
!     of time derivative operator: \partial / \partial t.
!
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix

  implicit none

!
! > Shift in time
! 
  ! fvEqn%xoo(:) = fvEqn%xo(:)
  ! fvEqn%yoo(:) = fvEqn%yo(:)
  ! fvEqn%zoo(:) = fvEqn%zo(:)

  ! fvEqn%xo(:) = U%x(:)
  ! fvEqn%yo(:) = U%y(:)
  ! fvEqn%zo(:) = U%z(:)

  ! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
  do inp=1,numCells

    !.....for u  sp => spu; for v  sp => spv; for w  sp => sp
    !.....sum source terms
    spu(inp) = 0.0_dp
    spv(inp) = 0.0_dp
    sp(inp)  = 0.0_dp

    !=======================================================================
    ! Unsteady term
    !=======================================================================
    if(bdf) then
    !-----------------------------------------------------------------------
    ! Backward differentiation formula:
    ! in case that BTIME=0. --> Implicit Euler
    ! in case that BTIME=1. --> Three Level Implicit Time Integration Method or BDF2
    !-----------------------------------------------------------------------
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

      spu(inp) = spu(inp) + apotime*(1+0.5*btime)
      spv(inp) = spv(inp) + apotime*(1+0.5*btime)
      sp(inp)  = sp(inp)  + apotime*(1+0.5*btime)
    !-----------------------------------------------------------------------
    endif

  end do


end subroutine


subroutine fvm_ddt_scalar_field
!  
!******************************************************************************
!
!     Finds source and matrix coefficients representing FVM discretization 
!     of time derivative operator: \partial / \partial t.
!
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix

  implicit none

!
! > Shift in time
! 
  ! fvEqn%oo(:) = fvEqn%o(:)
  ! fvEqn%o(:) = phi%mag(:)

  ! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
  do inp=1,numCells

    !.....for u  sp => spu; for v  sp => spv; for w  sp => sp
    !.....sum source terms
    spu(inp) = 0.0_dp

    !=======================================================================
    ! Unsteady term
    !=======================================================================
    if(bdf) then
    !-----------------------------------------------------------------------
    ! Backward differentiation formula:
    ! in case that BTIME=0. --> Implicit Euler
    ! in case that BTIME=1. --> Three Level Implicit Time Integration Method or BDF2
    !-----------------------------------------------------------------------
      apotime=den(inp)*vol(inp)/timestep

      sut = apotime*((1+btime)*uo(inp))
    
      if (btime > 0.99) then ! bdf2 scheme
        sut = sut - apotime*(0.5*btime*uoo(inp))
      endif

      su(inp) = su(inp) + sut

      spu(inp) = spu(inp) + apotime*(1+0.5*btime)
    !-----------------------------------------------------------------------
    endif

  end do


end subroutine



subroutine fvm_d2dt2_scalar_field
!  
!******************************************************************************
!
!   Description: 
!     Second order time differentiation for scalar fields.
!
!     Finds source and matrix coefficients representing FVM discretization 
!     of scnd. order time derivative operator: \partial^2 / \partial t^2.
!
!     System of linear equations is written as:
!     $ a_p^{(i)}*\phi_p^(i)-\sum_{j=1}^{nb} a_j^{(i)}*\phi_j^(i) = b_p{(i)}, i=1,ncells $
!
!******************************************************************************
!
  use types
  use parameters
  use geometry
  use sparse_matrix

  implicit none

!
! > Shift in time
! 
  ! fvEqn%oo(:) = fvEqn%o(:)
  ! fvEqn%o(:) = phi%mag(:)

  ! spu(1:numCells) = den(1:numCells)*vol(1:numCells)/timestep**2
  ! su(1:numCells) = den(1:numCells)*vol(1:numCells)/timestep**2 * (2*uo(1:numCells) - uoo(1:numCells))

  ! CALCULATE SOURCE TERMS INTEGRATED OVER VOLUME
  do inp=1,numCells

    !.....for u  sp => spu; for v  sp => spv; for w  sp => sp
    !.....sum source terms
    spu(inp) = 0.0_dp

    !=======================================================================
    ! Unsteady term
    !=======================================================================
    if(bdf) then
    !-----------------------------------------------------------------------
    ! Backward differentiation formula:
    ! in case that BTIME=0. --> Implicit Euler
    ! in case that BTIME=1. --> Three Level Implicit Time Integration Method or BDF2
    !-----------------------------------------------------------------------
      apotime=den(inp)*vol(inp)/timestep**2

      sut = apotime*(2*uo(inp) - uoo(inp))

      su(inp) = su(inp) + sut

      spu(inp) = spu(inp) + apotime
    !-----------------------------------------------------------------------
    endif

  end do


end subroutine