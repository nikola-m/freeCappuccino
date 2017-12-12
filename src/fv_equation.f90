module fv_equation
!
! Definition of Finite Volume Method linear equation object, and oparations over it.
!
use types
use geometry
use sparse_matrix
use tensor_fields

implicit none

!
! The fvEquation derived data type (all objects from fvm_discretisation_module belong to this type.)
!
type, extends(csrMatrix) :: fvEquation
  real(dp), dimension(:), allocatable :: source
  real(dp), dimension(:), allocatable :: o   ! Past time values of solution vector field 
  real(dp), dimension(:), allocatable :: oo  ! Second level past times for BDF2 time differencing scheme
end type fvEquation



type, extends(csrMatrix) :: fvVectorEquation

  ! Source terms that go into main diagonal
  real(dp), dimension(:), allocatable :: spu
  real(dp), dimension(:), allocatable :: spv
  real(dp), dimension(:), allocatable :: spw

  ! Sources
  real(dp), dimension(:), allocatable :: su
  real(dp), dimension(:), allocatable :: sv
  real(dp), dimension(:), allocatable :: sw

  ! Past time values of solution vector field 
  real(dp), dimension(:), allocatable :: xo
  real(dp), dimension(:), allocatable :: yo
  real(dp), dimension(:), allocatable :: zo


  real(dp), dimension(:), allocatable :: xoo
  real(dp), dimension(:), allocatable :: yoo
  real(dp), dimension(:), allocatable :: zoo

end type fvVectorEquation


interface operator(==)
   module procedure add_source_to_fvEquation
   module procedure add_volVectorFieldSource_to_fvVectorEquation
   module procedure add_fvEquations
   module procedure add_fvVectorEquations
end interface


! Overload summation to be able to add fvEquation and volScalarField
! This enables to have fvm_... routines which usually return FvEquation on rhs of == sign
interface operator(+)
   module procedure add_source_to_fvEquation  
   module procedure add_volVectorFieldSource_to_fvVectorEquation
   module procedure add_fvEquations
   module procedure add_fvVectorEquations
end interface


interface operator(-)
   module procedure substract_source_from_fvEquation  
   module procedure substract_volVectorFieldSource_from_fvVectorEquation
   module procedure substract_fvEquations
   module procedure substract_fvVectorEquations
end interface

public

contains



function new_csrMatrix( )
    implicit none
    integer :: i
    type(csrMatrix) :: new_csrMatrix
    
    allocate(new_csrMatrix%ioffset ( numCells+1 ))
    allocate(new_csrMatrix%ja ( nnz ))
    allocate(new_csrMatrix%coef ( nnz ))

    do i=1,numCells+1
      new_csrMatrix % ioffset(i) = ioffset(i)
    enddo
    do i=1,nnz
      new_csrMatrix % ja(i) = ja(i)
    enddo

end function new_csrMatrix


function new_fvEquation( ) result(fvEqn)
    implicit none
    integer :: i
    type(fvEquation) :: fvEqn

    allocate(fvEqn % ioffset ( numCells+1 ))
    allocate(fvEqn % ja ( nnz ))
    allocate(fvEqn % coef ( nnz ))

    allocate(fvEqn % source ( numCells ))

    allocate(fvEqn % o ( numTotal ))
    allocate(fvEqn % oo ( numTotal ))

    do i=1,numCells+1
      fvEqn % ioffset(i) = ioffset(i)
    enddo
    do i=1,nnz
      fvEqn % ja(i) = ja(i)
    enddo

end function new_fvEquation

function new_fvVectorEquation( ) result(fvEqn)
    implicit none
    integer :: i
    type(fvVectorEquation) :: fvEqn

    allocate(fvEqn % ioffset ( numCells+1 ))
    allocate(fvEqn % ja ( nnz ))
    allocate(fvEqn % coef ( nnz ))

    allocate(fvEqn % su ( numCells ))
    allocate(fvEqn % sv ( numCells ))
    allocate(fvEqn % sw ( numCells ))

    allocate(fvEqn % spu ( numCells ))
    allocate(fvEqn % spv ( numCells ))
    allocate(fvEqn % spw ( numCells ))

    allocate(fvEqn % xo ( numTotal ))
    allocate(fvEqn % yo ( numTotal ))
    allocate(fvEqn % zo ( numTotal ))

    allocate(fvEqn % xoo ( numTotal ))
    allocate(fvEqn % yoo ( numTotal ))
    allocate(fvEqn % zoo ( numTotal ))

    do i=1,numCells+1
      fvEqn % ioffset(i) = ioffset(i)
    enddo
    do i=1,nnz
      fvEqn % ja(i) = ja(i)
    enddo

end function new_fvVectorEquation




function add_source_to_fvEquation(fvEqnIn,source) result( fvEqnOut )
!
! Adds source to eqn. system rhs vector.
!
  use types
  implicit none
!
! > Input
!
  type(fvEquation), intent(in) :: fvEqnIn
  type(volScalarField), intent(in) :: source
!
! > Result
!
  type(fvEquation) :: fvEqnOut

  integer :: i

  fvEqnOut = new_fvEquation()

  do i=1,numCells
      fvEqnOut % source(i) = fvEqnIn % source(i) + source % mag(i)
  enddo

end function add_source_to_fvEquation




function add_volVectorFieldSource_to_fvVectorEquation(fvEqnIn,vecSource) result( fvEqnOut )
!
! Adds source to eqn. system rhs vector.
!
  use types
  implicit none
!
! > Input
!
  type(fvVectorEquation), intent(in) :: fvEqnIn
  type(volVectorField), intent(in) :: vecSource
!
! > Result
!
  type(fvVectorEquation) :: fvEqnOut

  integer :: i

  fvEqnOut = new_fvVectorEquation()

  do i=1,numCells
      fvEqnOut % su(i) = fvEqnIn % su(i) + vecSource % x(i)
      fvEqnOut % sv(i) = fvEqnIn % sv(i) + vecSource % y(i)
      fvEqnOut % sw(i) = fvEqnIn % sw(i) + vecSource % z(i)
  enddo

end function add_volVectorFieldSource_to_fvVectorEquation



function substract_source_from_fvEquation(fvEqnIn,source) result( fvEqnOut )
!
! Adds source to eqn. system rhs vector.
!
  use types
  implicit none
!
! > Input
!
  type(fvEquation), intent(in) :: fvEqnIn
  type(volScalarField), intent(in) :: source
!
! > Result
!
  type(fvEquation) :: fvEqnOut

  integer :: i

  fvEqnOut = new_fvEquation()

  do i=1,numCells
      fvEqnOut % source(i) = fvEqnIn % source(i) - source % mag(i)
  enddo

end function substract_source_from_fvEquation




function substract_volVectorFieldSource_from_fvVectorEquation(fvEqn,vecSource) result( fvEqn_new )
!
! Adds source to eqn. system rhs vector.
!
  use types
  implicit none
!
! > Input
!
  type(fvVectorEquation), intent(in) :: fvEqn
  type(volVectorField), intent(in) :: vecSource
!
! > Result
!
  type(fvVectorEquation) :: fvEqn_new

  integer :: i

  fvEqn_new = new_fvVectorEquation()

  do i=1,numCells
      fvEqn_new % su(i) = fvEqn % su(i) - vecSource % x(i)
      fvEqn_new % sv(i) = fvEqn % sv(i) - vecSource % y(i)
      fvEqn_new % sw(i) = fvEqn % sw(i) - vecSource % z(i)
  enddo

end function substract_volVectorFieldSource_from_fvVectorEquation


function add_fvEquations(fvEqn1,fvEqn2) result( fvEqn_new )
!
! Adds two objects of type(fvEquation).
!
  use types
  implicit none
!
! > Input
!
  type(fvEquation), intent(in) :: fvEqn1
  type(fvEquation), intent(in) :: fvEqn2
!
! > Result
!
  type(fvEquation) :: fvEqn_new

  integer :: i

  fvEqn_new = new_fvEquation()

  do i=1,numCells
      fvEqn_new % source(i) = fvEqn1 % source(i) + fvEqn2 % source(i)
  enddo
  do i=1,nnz
      fvEqn_new % coef(i) = fvEqn1 % coef(i) + fvEqn2 % coef(i)
  enddo

end function add_fvEquations


function add_fvVectorEquations(fvEqn1,fvEqn2) result( fvEqn_new )
!
! Adds two objects of type(fvVectorEquation).
!
  use types
  implicit none
!
! > Input
!
  type(fvVectorEquation), intent(in) :: fvEqn1
  type(fvVectorEquation), intent(in) :: fvEqn2
!
! > Result
!
  type(fvVectorEquation) :: fvEqn_new

  integer :: i

  fvEqn_new = new_fvVectorEquation()

  do i=1,numCells
      fvEqn_new % su(i) = fvEqn1 % su(i) + fvEqn2 % su(i)
      fvEqn_new % sv(i) = fvEqn1 % sv(i) + fvEqn2 % sv(i)
      fvEqn_new % sw(i) = fvEqn1 % sw(i) + fvEqn2 % sw(i)

      fvEqn_new % spu(i) = fvEqn1 % spu(i) + fvEqn2 % spu(i)
      fvEqn_new % spv(i) = fvEqn1 % spv(i) + fvEqn2 % spv(i)
      fvEqn_new % spw(i) = fvEqn1 % spw(i) + fvEqn2 % spw(i)
  enddo
  do i=1,nnz
      fvEqn_new % coef(i) = fvEqn1 % coef(i) + fvEqn2 % coef(i)
  enddo

end function add_fvVectorEquations


function substract_fvEquations(fvEqn1,fvEqn2) result( fvEqn_new )
!
! Substracts two objects of type(fvEquation).
!
  use types
  implicit none
!
! > Input
!
  type(fvEquation), intent(in) :: fvEqn1
  type(fvEquation), intent(in) :: fvEqn2
!
! > Result
!
  type(fvEquation) :: fvEqn_new

  integer :: i

  fvEqn_new = new_fvEquation()

  do i=1,numCells
      fvEqn_new % source(i) = fvEqn1 % source(i) - fvEqn2 % source(i)
  enddo
  do i=1,nnz
      fvEqn_new % coef(i) = fvEqn1 % coef(i) - fvEqn2 % coef(i)
  enddo

end function substract_fvEquations


function substract_fvVectorEquations(fvEqn1,fvEqn2) result( fvEqn )
!
! Substracts two objects of type(fvVectorEquation).
!
  use types
  implicit none
!
! > Input
!
  type(fvVectorEquation), intent(in) :: fvEqn1
  type(fvVectorEquation), intent(in) :: fvEqn2
!
! > Result
!
  type(fvVectorEquation) :: fvEqn

  integer :: i

  fvEqn = new_fvVectorEquation()

  do i=1,numCells
      fvEqn % su(i) = fvEqn1 % su(i) - fvEqn2 % su(i)
      fvEqn % sv(i) = fvEqn1 % sv(i) - fvEqn2 % sv(i)
      fvEqn % sw(i) = fvEqn1 % sw(i) - fvEqn2 % sw(i)

      fvEqn % spu(i) = fvEqn1 % spu(i) - fvEqn2 % spu(i)
      fvEqn % spv(i) = fvEqn1 % spv(i) - fvEqn2 % spv(i)
      fvEqn % spw(i) = fvEqn1 % spw(i) - fvEqn2 % spw(i)
  enddo
  do i=1,nnz
      fvEqn % coef(i) = fvEqn1 % coef(i) - fvEqn2 % coef(i)
  enddo

end function substract_fvVectorEquations


end module fv_equation
