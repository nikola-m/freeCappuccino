module tensor_fields
!
! Definition of volume and surface tensor fields and operations on them.
!
use types
use geometry, only: numCells,numInnerFaces
implicit none


!
! > Volume fields
!

type volScalarField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: mag
end type

type volVectorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: x, y, z
end type

type volSymmetricTensorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: xx, xy, xz
  real(dp), dimension(:), allocatable ::     yy, yz
  real(dp), dimension(:), allocatable ::         zz
end type

type volTensorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: xx, xy, xz
  real(dp), dimension(:), allocatable :: yx, yy, yz
  real(dp), dimension(:), allocatable :: zx, zy, zz
end type


!
! > Surface fields
!

type surfaceScalarField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: mag
end type

type surfaceVectorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: x, y, z
end type

type surfaceSymmetricTensorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: xx, xy, xz
  real(dp), dimension(:), allocatable ::     yy, yz
  real(dp), dimension(:), allocatable ::         zz
end type

type surfaceTensorField
  character(len=30) :: field_name
  real(dp), dimension(:), allocatable :: xx, xy, xz
  real(dp), dimension(:), allocatable :: yx, yy, yz
  real(dp), dimension(:), allocatable :: zx, zy, zz
end type


!
! > Operations on vector and tensor fields
!

interface operator(.dot.)
   module procedure calc_inner_product
   module procedure calc_inner_product_surface_vectors
   module procedure calc_inner_product_rank2_and_rank1_tensors
   module procedure calc_inner_product_rank1_and_rank2_tensors
end interface

interface operator(.cross.)
   module procedure calc_cross_product
end interface

interface operator(.tensor.)
   module procedure calc_tensor_product
   module procedure calc_tensor_product_rank2_tensors
end interface

interface operator(.colon.)
   module procedure calc_inner_product_rank2_tensor
   module procedure calc_inner_product_rank2_symmetric_tensor
end interface

interface operator(.transposed.)
   module procedure transpose_rank2_tensor
end interface

interface operator(.trace.)
   module procedure trace_rank2_tensor
   module procedure trace_rank2_symmetric_tensor
end interface

interface operator(.det.)
   module procedure determinant_rank2_tensor
   module procedure determinant_rank2_symmetric_tensor
end interface

interface operator(.diagonal.)
   module procedure diagonal
end interface

interface operator(.hodge.)
   module procedure hodge_dual
end interface

interface operator(.curl.)
   module procedure curl
end interface

interface operator(.symm.)
   module procedure symm
end interface

interface operator(.skew.)
   module procedure skew
end interface

interface operator(.magSq.)
   module procedure magSqrTensorField
   module procedure magSqrSymmetricTensorField
end interface

interface operator(.mag.)
   module procedure magTensorField
   module procedure magSymmetricTensorField
end interface

interface operator(.dev.)
   module procedure deviatoric_part_rank2_tensor
end interface

interface operator(.devTwo.)
   module procedure deviatoric_part_rank2_tensor_23
end interface

interface operator(.hyd.)
   module procedure hydrostatic_part_rank2_tensor
end interface



! > Operator overloading for vector and tensor fields

interface operator(+)
   module procedure add_tensors
end interface

interface operator(-)
   module procedure substract_tensors
end interface

interface operator(*)
   module procedure scalar_vector_multiply
   module procedure scalar_field_vector_multiply
   module procedure scalar_rank2_tensor_multiply
   module procedure scalar_field_rank2_tensor_multiply
   module procedure scalar_surface_vector_multiply
   module procedure scalar_field_surface_vector_multiply
end interface


public


contains


!
! > Create new fields
!

function new_volScalarField( num )
    implicit none
    integer, intent(in) :: num
    type(volScalarField) :: new_volScalarField

    new_volScalarField % field_name = 'unknownVolScalarField'

    allocate(new_volScalarField%mag ( num ))
end function new_volScalarField


function new_volVectorField( num )
    implicit none
    integer, intent(in) :: num
    type(volVectorField) :: new_volVectorField

    new_volVectorField % field_name = 'unknownVolVectorField'

    allocate(new_volVectorField%x ( num ))
    allocate(new_volVectorField%y ( num ))
    allocate(new_volVectorField%z ( num ))
end function new_volVectorField


function new_volSymmetricTensorField( num )
    implicit none
    integer, intent(in) :: num
    type(volSymmetricTensorField) :: new_volSymmetricTensorField

    new_volSymmetricTensorField % field_name = 'unknownVolSymmTensorField'

    allocate(new_volSymmetricTensorField%xx ( num ))
    allocate(new_volSymmetricTensorField%xy ( num ))
    allocate(new_volSymmetricTensorField%xz ( num ))

    allocate(new_volSymmetricTensorField%yy ( num ))
    allocate(new_volSymmetricTensorField%yz ( num ))

    allocate(new_volSymmetricTensorField%zz ( num ))
end function new_volSymmetricTensorField

function new_volTensorField( num )
    implicit none
    integer, intent(in) :: num
    type(volTensorField) :: new_volTensorField

    new_volTensorField % field_name = 'unknownVolTensorField'

    allocate(new_volTensorField%xx ( num ))
    allocate(new_volTensorField%xy ( num ))
    allocate(new_volTensorField%xz ( num ))

    allocate(new_volTensorField%yx ( num ))
    allocate(new_volTensorField%yy ( num ))
    allocate(new_volTensorField%yz ( num ))

    allocate(new_volTensorField%zx ( num ))
    allocate(new_volTensorField%zy ( num ))
    allocate(new_volTensorField%zz ( num ))
end function new_volTensorField

!

function new_surfaceScalarField( num )
    implicit none
    integer, intent(in) :: num
    type(surfaceScalarField) :: new_surfaceScalarField

    new_surfaceScalarField % field_name = 'unknownSurfaceScalarField'

    allocate(new_surfaceScalarField%mag ( num ))
end function new_surfaceScalarField

function new_surfaceVectorField(num)
    implicit none
    integer, intent(in) :: num
    type(surfaceVectorField) :: new_surfaceVectorField

    new_surfaceVectorField % field_name = 'unknownSurfaceVectorField'

    allocate(new_surfaceVectorField%x ( num ))
    allocate(new_surfaceVectorField%y ( num ))
    allocate(new_surfaceVectorField%z ( num ))
end function new_surfaceVectorField

function new_surfaceSymmetricTensorField( num )
    implicit none
    integer, intent(in) :: num
    type(surfaceSymmetricTensorField) :: new_surfaceSymmetricTensorField

    new_surfaceSymmetricTensorField % field_name = 'unknownSurfSymTensorField'

    allocate(new_surfaceSymmetricTensorField%xx ( num ))
    allocate(new_surfaceSymmetricTensorField%xy ( num ))
    allocate(new_surfaceSymmetricTensorField%xz ( num ))

    allocate(new_surfaceSymmetricTensorField%yy ( num ))
    allocate(new_surfaceSymmetricTensorField%yz ( num ))

    allocate(new_surfaceSymmetricTensorField%zz ( num ))
end function new_surfaceSymmetricTensorField

function new_surfaceTensorField(num)
    implicit none
    integer, intent(in) :: num
    type(surfaceTensorField) :: new_surfaceTensorField

    new_surfaceTensorField % field_name = 'unknownSurfaceTensorField'

    allocate(new_surfaceTensorField%xx ( num ))
    allocate(new_surfaceTensorField%xy ( num ))
    allocate(new_surfaceTensorField%xz ( num ))

    allocate(new_surfaceTensorField%yx ( num ))
    allocate(new_surfaceTensorField%yy ( num ))
    allocate(new_surfaceTensorField%yz ( num ))

    allocate(new_surfaceTensorField%xz ( num ))
    allocate(new_surfaceTensorField%yz ( num ))
    allocate(new_surfaceTensorField%zz ( num ))
end function new_surfaceTensorField

!
! > Operations over fields
!

! The .dot. operator defining scalar product between two vector fields and ...

function calc_inner_product(v1, v2)  result(inner_product)
    implicit none
    type(volVectorField), intent(in) :: v1, v2
    real(dp), dimension(numCells) :: inner_product
    integer :: i

    do i = 1,numCells
        inner_product(i) = v1%x(i) * v2%x(i) + v1%y(i) * v2%y(i) + v1%z(i) * v2%z(i)
    enddo
end function calc_inner_product

function calc_inner_product_surface_vectors(v1, v2)  result(inner_product)
    implicit none
    type(surfaceVectorField), intent(in) :: v1, v2
    real(dp), dimension(numInnerFaces) :: inner_product
    integer :: i

    do i = 1,numInnerFaces
        inner_product(i) = v1%x(i) * v2%x(i) + v1%y(i) * v2%y(i) + v1%z(i) * v2%z(i)
    enddo
end function calc_inner_product_surface_vectors

! ... inner product between tensor and vector vi = Tij*vj, and...

function calc_inner_product_rank2_and_rank1_tensors(T1, v1)  result(v2)
    implicit none
    type(volTensorField), intent(in) :: T1
    type(volVectorField), intent(in) :: v1
    type(volVectorField)             :: v2
    integer :: i

    v2 = new_volVectorField(numCells)

    do i = 1,numCells
        v2%x (i) = T1%xx(i) * v1%x(i) + T1%xy(i) * v1%y(i) + T1%xz(i) * v1%z(i)  
        v2%y (i) = T1%yx(i) * v1%x(i) + T1%yy(i) * v1%y(i) + T1%yz(i) * v1%z(i)  
        v2%z (i) = T1%zx(i) * v1%x(i) + T1%zy(i) * v1%y(i) + T1%zz(i) * v1%z(i)
    enddo
end function calc_inner_product_rank2_and_rank1_tensors

! ... inner product between vector and tensor vi = vj*Tij

function calc_inner_product_rank1_and_rank2_tensors(v1,T1)  result(v2)
    implicit none
    type(volVectorField), intent(in) :: v1
    type(volTensorField), intent(in) :: T1
    type(volVectorField)             :: v2

    v2 = new_volVectorField(numCells)
!
!   bi = Tji*aj using derived operators:
!
!        transpose tensor Tij
!        |               inner product bi = Tij*aj
!        |               |
    v2 = .transposed.T1 .dot. v1

end function calc_inner_product_rank1_and_rank2_tensors


! The cross product of two vector fields

function calc_cross_product(v1, v2)  result(v3)
    implicit none
    type(volVectorField), intent(in) :: v1, v2
    type(volVectorField)             :: v3
    integer :: i

    v3 = new_volVectorField(numCells)

    do i = 1,numCells
        v3%x (i) = v1%y(i) * v2%z(i) - v1%z(i) * v2%y(i)  
        v3%y (i) = v1%z(i) * v2%x(i) - v1%x(i) * v2%z(i) 
        v3%z (i) = v1%x(i) * v2%y(i) - v1%y(i) * v2%x(i)
    enddo
end function calc_cross_product


! ! The .tensor. operator defining tensor, or outer product between two column vectors, or ...

function calc_tensor_product(v1, v2)  result(T1)
    implicit none
    type(volVectorField), intent(in) :: v1, v2
    type(volTensorField)             :: T1
    integer :: i

    T1 = new_volTensorField(numCells)

    do i = 1,numCells
        T1%xx(i) = v1%x(i) * v2%x(i) 
        T1%xy(i) = v1%x(i) * v2%y(i) 
        T1%xz(i) = v1%x(i) * v2%z(i) 

        T1%yx(i) = v1%y(i) * v2%x(i) 
        T1%yy(i) = v1%y(i) * v2%y(i) 
        T1%yz(i) = v1%y(i) * v2%z(i)

        T1%zx(i) = v1%z(i) * v2%x(i) 
        T1%zy(i) = v1%z(i) * v2%y(i) 
        T1%zz(i) = v1%z(i) * v2%z(i)
    enddo
end function calc_tensor_product

! ... between two tensors

function calc_tensor_product_rank2_tensors(T1, T2)  result(T3)
    implicit none
    type(volTensorField), intent(in) :: T1, T2
    type(volTensorField)             :: T3
    integer :: i

    T3 = new_volTensorField(numCells)

    do i = 1,numCells
        T3%xx(i) = T1%xx(i) * T2%xx(i) + T1%xy(i) * T2%yx(i) + T1%xz(i) * T2%zx(i)
        T3%xy(i) = T1%xx(i) * T2%xy(i) + T1%xy(i) * T2%yy(i) + T1%xz(i) * T2%zy(i) 
        T3%xz(i) = T1%xx(i) * T2%xz(i) + T1%xy(i) * T2%yz(i) + T1%xz(i) * T2%zz(i) 

        T3%yx(i) = T1%yx(i) * T2%xx(i) + T1%yy(i) * T2%yx(i) + T1%yz(i) * T2%zx(i) 
        T3%yy(i) = T1%yx(i) * T2%xy(i) + T1%yy(i) * T2%yy(i) + T1%yz(i) * T2%zy(i) 
        T3%yz(i) = T1%yx(i) * T2%xz(i) + T1%yy(i) * T2%yz(i) + T1%yz(i) * T2%zz(i)

        T3%zx(i) = T1%zx(i) * T2%xx(i) + T1%zx(i) * T2%yx(i) + T1%zz(i) * T2%zx(i)
        T3%zy(i) = T1%zx(i) * T2%xy(i) + T1%zy(i) * T2%yy(i) + T1%zz(i) * T2%zy(i) 
        T3%zz(i) = T1%zx(i) * T2%xz(i) + T1%zy(i) * T2%yz(i) + T1%zz(i) * T2%zz(i)
    enddo
end function calc_tensor_product_rank2_tensors



function calc_inner_product_rank2_tensor(T1, T2) result(inner_product)
    implicit none
    type(volTensorField), intent(in) :: T1, T2
    type(volScalarField) :: inner_product
    integer :: i

    inner_product = new_volScalarField(numCells)

    do i = 1,numCells
        inner_product%mag(i) = T1%xx(i) * T2%xx(i) + T1%xy(i) * T2%xy(i) + T1%xz(i) * T2%xz(i)  &
                             + T1%yx(i) * T2%yx(i) + T1%yy(i) * T2%yy(i) + T1%yz(i) * T2%yz(i)  &
                             + T1%zx(i) * T2%zx(i) + T1%zy(i) * T2%zy(i) + T1%zz(i) * T2%zz(i)
    enddo
end function calc_inner_product_rank2_tensor

function calc_inner_product_rank2_symmetric_tensor(T1, T2) result(inner_product)
    implicit none
    type(volSymmetricTensorField), intent(in) :: T1, T2
    type(volScalarField) :: inner_product
    integer :: i

    inner_product = new_volScalarField(numCells)

    do i = 1,numCells
        inner_product%mag(i) = T1%xx(i) * T2%xx(i) + T1%xy(i) * T2%xy(i) + T1%xz(i) * T2%xz(i)  &
                                                   + T1%yy(i) * T2%yy(i) + T1%yz(i) * T2%yz(i)  &
                                                                         + T1%zz(i) * T2%zz(i)
    enddo
end function calc_inner_product_rank2_symmetric_tensor


function transpose_rank2_tensor(T1)  result(T2)
    implicit none
    type(volTensorField), intent(in) :: T1
    type(volTensorField)             :: T2
    integer :: i

    T2 = new_volTensorField(numCells)

    do i = 1,numCells
        T2%xx(i) = T1%xx(i)
        T2%xy(i) = T1%yx(i) 
        T2%xz(i) = T1%zx(i)

        T2%yx(i) = T1%xy(i)
        T2%yy(i) = T1%yy(i)
        T2%yz(i) = T1%zy(i)

        T2%zx(i) = T1%xz(i)
        T2%zy(i) = T1%yz(i)  
        T2%zz(i) = T1%zz(i)
    enddo
end function transpose_rank2_tensor

function add_tensors(T1,T2)  result(T3)
    implicit none
    type(volTensorField), intent(in) :: T1, T2
    type(volTensorField)             :: T3
    integer :: i

    T3 = new_volTensorField(numCells)

    do i = 1,numCells
        T3%xx(i) = T1%xx(i) + T2%xx(i)
        T3%xy(i) = T1%xy(i) + T2%xy(i)
        T3%xz(i) = T1%xz(i) + T2%xz(i)

        T3%yx(i) = T1%yx(i) + T2%yx(i)
        T3%yy(i) = T1%yy(i) + T2%yy(i)
        T3%yz(i) = T1%yz(i) + T2%yz(i)

        T3%zx(i) = T1%zx(i) + T2%zx(i)
        T3%zy(i) = T1%zy(i) + T2%zy(i)
        T3%zz(i) = T1%zz(i) + T2%zz(i)
    enddo
end function add_tensors

function substract_tensors(T1,T2)  result(T3)
    implicit none
    type(volTensorField), intent(in) :: T1, T2
    type(volTensorField)             :: T3
    integer :: i

    T3 = new_volTensorField(numCells)

    do i = 1,numCells
        T3%xx(i) = T1%xx(i) - T2%xx(i)
        T3%xy(i) = T1%xy(i) - T2%xy(i)
        T3%xz(i) = T1%xz(i) - T2%xz(i)

        T3%yx(i) = T1%yx(i) - T2%yx(i)
        T3%yy(i) = T1%yy(i) - T2%yy(i)
        T3%yz(i) = T1%yz(i) - T2%yz(i)

        T3%zx(i) = T1%zx(i) - T2%zx(i)
        T3%zy(i) = T1%zy(i) - T2%zy(i)
        T3%zz(i) = T1%zz(i) - T2%zz(i)
    enddo
end function substract_tensors


function scalar_vector_multiply(alpha,v1)  result(v2)
    implicit none
    real(dp), intent(in) :: alpha
    type(volVectorField), intent(in) :: v1
    type(volVectorField)             :: v2
    integer :: i

    v2 = new_volVectorField(numCells)

    do i = 1,numCells
        v2%x (i) = alpha * v1%x(i)  
        v2%y (i) = alpha * v1%y(i) 
        v2%z (i) = alpha * v1%z(i)
    enddo
end function scalar_vector_multiply


function scalar_field_vector_multiply(alpha,v1)  result(v2)
    implicit none
    real(dp), dimension(numCells), intent(in) :: alpha
    type(volVectorField), intent(in) :: v1
    type(volVectorField)             :: v2
    integer :: i

    v2 = new_volVectorField(numCells)

    do i = 1,numCells
        v2%x (i) = alpha(i) * v1%x(i)  
        v2%y (i) = alpha(i) * v1%y(i) 
        v2%z (i) = alpha(i) * v1%z(i)
    enddo
end function scalar_field_vector_multiply


function scalar_rank2_tensor_multiply(alpha,T1)  result(T2)
    implicit none
    real(dp), intent(in) :: alpha
    type(volTensorField), intent(in) :: T1
    type(volTensorField)             :: T2
    integer :: i

    T2 = new_volTensorField(numCells)

    do i = 1,numCells
        T2%xx(i) = alpha * T1%xx(i)
        T2%xy(i) = alpha * T1%xy(i)
        T2%xz(i) = alpha * T1%xz(i)

        T2%yx(i) = alpha * T1%yx(i)
        T2%yy(i) = alpha * T1%yy(i)
        T2%yz(i) = alpha * T1%yz(i)

        T2%zx(i) = alpha * T1%zx(i)
        T2%zy(i) = alpha * T1%zy(i)
        T2%zz(i) = alpha * T1%zz(i)
    enddo
end function scalar_rank2_tensor_multiply


function scalar_field_rank2_tensor_multiply(alpha,T1)  result(T2)
    implicit none
    real(dp), dimension(numCells), intent(in) :: alpha
    type(volTensorField), intent(in) :: T1
    type(volTensorField)             :: T2
    integer :: i

    T2 = new_volTensorField(numCells)

    do i = 1,numCells
        T2%xx(i) = alpha(i) * T1%xx(i)
        T2%xy(i) = alpha(i) * T1%xy(i)
        T2%xz(i) = alpha(i) * T1%xz(i)

        T2%yx(i) = alpha(i) * T1%yx(i)
        T2%yy(i) = alpha(i) * T1%yy(i)
        T2%yz(i) = alpha(i) * T1%yz(i)

        T2%zx(i) = alpha(i) * T1%zx(i)
        T2%zy(i) = alpha(i) * T1%zy(i)
        T2%zz(i) = alpha(i) * T1%zz(i)
    enddo
end function scalar_field_rank2_tensor_multiply


function scalar_surface_vector_multiply(alpha,v1)  result(v2)
    implicit none
    real(dp), intent(in) :: alpha
    type(surfaceVectorField), intent(in) :: v1
    type(surfaceVectorField)             :: v2
    integer :: i

    v2 = new_surfaceVectorField(numCells)

    do i = 1,numCells
        v2%x (i) = alpha * v1%x(i)  
        v2%y (i) = alpha * v1%y(i) 
        v2%z (i) = alpha * v1%z(i)
    enddo
end function scalar_surface_vector_multiply


function scalar_field_surface_vector_multiply(alpha,v1)  result(v2)
    implicit none
    real(dp), dimension(numCells), intent(in) :: alpha
    type(surfaceVectorField), intent(in) :: v1
    type(surfaceVectorField)             :: v2
    integer :: i

    v2 = new_surfaceVectorField(numCells)

    do i = 1,numCells
        v2%x (i) = alpha(i) * v1%x(i)  
        v2%y (i) = alpha(i) * v1%y(i) 
        v2%z (i) = alpha(i) * v1%z(i)
    enddo
end function scalar_field_surface_vector_multiply


function trace_rank2_symmetric_tensor(T) result(trace)
    implicit none
    type(volSymmetricTensorField), intent(in) :: T
    real(dp), dimension(numCells) :: trace
    integer :: i

    do i = 1,numCells
        trace(i) = T%xx(i) + T%yy(i) + T%zz(i)
    enddo
end function trace_rank2_symmetric_tensor


function trace_rank2_tensor(T) result(trace)
    implicit none
    type(volTensorField), intent(in) :: T
    real(dp), dimension(numCells) :: trace
    integer :: i

    do i = 1,numCells
        trace(i) = T%xx(i) + T%yy(i) + T%zz(i)
    enddo
end function trace_rank2_tensor


function determinant_rank2_symmetric_tensor(T) result(determinant)
    implicit none
    type(volSymmetricTensorField), intent(in) :: T
    real(dp), dimension(numCells) :: determinant
    integer :: i

    do i = 1,numCells
        determinant(i) = ( T%xx(i) * ( T%yy(i) * T%zz(i) - T%yz(i)*T%yz(i) ) - &
                           T%xy(i) * ( T%xy(i) * T%zz(i) - T%yz(i)*T%xz(i) ) + &
                           T%xz(i) * ( T%xy(i) * T%yz(i) - T%yy(i)*T%xz(i) )   )
    enddo
end function determinant_rank2_symmetric_tensor


function determinant_rank2_tensor(T) result(determinant)
    implicit none
    type(volTensorField), intent(in) :: T
    real(dp), dimension(numCells) :: determinant
    integer :: i

    do i = 1,numCells
        determinant(i) = ( T%xx(i) * ( T%yy(i) * T%zz(i) - T%yz(i)*T%zy(i) ) - &
                           T%xy(i) * ( T%yx(i) * T%zz(i) - T%yz(i)*T%zx(i) ) + &
                           T%xz(i) * ( T%yx(i) * T%zy(i) - T%yy(i)*T%zx(i) )   )
    enddo
end function determinant_rank2_tensor


function diagonal(T) result(v)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volVectorField)             :: v
    integer :: i

    v = new_volVectorField(numCells)

    do i = 1,numCells
        v%x (i) = T%xx(i)  
        v%y (i) = T%yy(i) 
        v%z (i) = T%zz(i)
    enddo
end function diagonal


function hodge_dual(T) result(v)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volVectorField)             :: v
    integer :: i

    v = new_volVectorField(numCells)

    do i = 1,numCells
        v%x (i) = T%yz(i)  
        v%y (i) =-T%xz(i) 
        v%z (i) = T%xy(i)
    enddo
end function hodge_dual


function eye(num) result(T)
    implicit none
    integer, intent(in) :: num
    type(volTensorField) :: T
    integer :: i

    T = new_volTensorField(num)

    do i = 1,num
        T%xx(i) = 1.0_dp  
        T%xy(i) = 0.0_dp  
        T%xz(i) = 0.0_dp  

        T%yx(i) = 0.0_dp 
        T%yy(i) = 1.0_dp 
        T%yz(i) = 0.0_dp 
 
        T%zx(i) = 0.0_dp 
        T%zy(i) = 0.0_dp 
        T%zz(i) = 1.0_dp
    enddo
end function eye

function symm(T)  result(D)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volTensorField)             :: D

    D = new_volTensorField(numCells)
!              overloaded operator - here '*' multiplies tensor fields by a constant scalar        
!              |     overloaded operator - here '+' adds two tensor fields
!              |     |
    D = 0.5_dp * ( T + .transposed.T )

end function symm

function skew(T)  result(S)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volTensorField)             :: S

    !S = new_volTensorField(numCells)
!              overloaded operator - here '*' multiplies tensor field by a constant scalar        
!              |     overloaded operator - here '-' substracts two tensor fields
!              |     |
    S = 0.5_dp * ( T - .transposed.T )

end function skew

function curl(D) result(v)
!
! Curl of a vector field is twice Hodge dual of skew-symmetric part of gradient tensor of that vector field.
! Make sure that D in this function call is gradient tensor, e.g. velocity gradient tensor, so it makes sense.
! To obtain it you will have to use gradient function from teh 'fvc' module, and apply it to a vector field
! in question: [volTensorField] D = fvc_grad([volVectorField] u)
!
    implicit none
    type(volTensorField), intent(in) :: D
    type(volVectorField)             :: v

    ! Will not allocate it because it will be allocated in .hodge. funtion
    !v = new_volVectorField(numCells)
!              overloaded operator - here '*' multiplies vector field by a constant scalar              
!              |  derived operators - Hodge dual of skew-symmetric part of tensor T
!              |  |            
    v = 2.0_dp * (.hodge.(.skew.D))
    
end function curl


function deviatoric_part_rank2_tensor(T)  result(devT)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volTensorField)             :: devT
    type(volTensorField)             :: I

    devT = new_volTensorField(numCells)
    I  = eye(numCells)
!            overloaded operator - here '-' substracts two tensor fields
!            |             overloaded operator - here '*' multiplies tensor field by a constant scalar
!            |             |            overloaded operator - here '*' multiplies tensor fields by a real scalar array of size[1:numCells]
!            |             |            |
    devT = T - ( 1./3.0_dp * ( .trace.T * I) )

end function deviatoric_part_rank2_tensor


function deviatoric_part_rank2_tensor_23(T)  result(devT)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volTensorField)             :: devT
    type(volTensorField)             :: I

    devT = new_volTensorField(numCells)
    I  = eye(numCells)
!            overloaded operator - here '-' substracts two tensor fields
!            |             overloaded operator - here '*' multiplies tensor field by a constant scalar
!            |             |            overloaded operator - here '*' multiplies tensor fields by a real scalar array of size[1:numCells]
!            |             |            |
    devT = T - ( 2./3.0_dp * ( .trace.T * I) )
                !^
                !!----2/3 here !

end function deviatoric_part_rank2_tensor_23


function hydrostatic_part_rank2_tensor(T)  result(hydT)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volTensorField)             :: hydT
    type(volTensorField)             :: I

    hydT = new_volTensorField(numCells)
    I  = eye(numCells)
!                     overloaded operator - here '*' multiplies tensor field by a constant scalar
!                     |            overloaded operator - here '*' multiplies tensor fields by a real scalar array of size[1:numCells]
!                     |            |
    hydT =  1./3.0_dp * ( .trace.T * I) 

end function hydrostatic_part_rank2_tensor

! .magSq.

function magSqrTensorField(T) result(scalar)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volScalarField)             :: scalar

        scalar = T.colon.T

end function

function magSqrSymmetricTensorField(S) result(scalar)
    implicit none
    type(volSymmetricTensorField), intent(in) :: S
    type(volScalarField)             :: scalar

        scalar = S.colon.S

end function

! .mag.

function magTensorField(T) result(scalar)
    implicit none
    type(volTensorField), intent(in) :: T
    type(volScalarField)             :: scalar

        scalar = T.colon.T
        scalar%mag = sqrt(scalar%mag)

end function

function magSymmetricTensorField(S) result(scalar)
    implicit none
    type(volSymmetricTensorField), intent(in) :: S
    type(volScalarField)             :: scalar

        scalar = S.colon.S
        scalar%mag = sqrt(scalar%mag)

end function

end module tensor_fields