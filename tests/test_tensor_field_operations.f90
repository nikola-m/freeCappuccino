  program testFieldOperation
  use types
  use utils
  use geometry
  use tensor_fields
  use output

  implicit none

! Locals
  integer :: ierr

  type(volVectorField) :: vec1, vec2
  type(volScalarField) :: phi
  type(volVectorField) :: resVec
  type(volTensorField) :: T



  ! 1) Print code logo and timestamp in monitor file
!+-----------------------------------------------------------------------------+

  call show_logo


  ! 2) Open & Read mesh file, calculate mesh geometrical quantities, allocate arrays
!+-----------------------------------------------------------------------------+

  call mesh_geometry


  ! 3) Define sparsity pattern according to given mesh connectivity data
!+-----------------------------------------------------------------------------+

  !call create_CSR_matrix_from_mesh_data


  ! 4) Finite Volume Discretization
!+-----------------------------------------------------------------------------+

  !allocate( phi(numCells) )
  ! vec1 = new_volVectorField( numCells )
  ! vec2 = new_volVectorField( numCells )
  phi = new_volScalarField( numCells )

!
! > Initialize vector fields
!

  vec1 = volVectorField(             &
                        "Vector1",   &
                        xc**2+yc+zc, &
                        xc+yc**2+zc, &
                        xc+yc+zc**2  &
                       )

  vec2 = volVectorField(           &
                        "Vector2", &
                        xc+yc+zc,  &
                        xc+yc+zc,  &
                        xc+yc+zc   &
                       )

!
! 1) Test dot product of two vector field
!
  phi%mag = vec1.dot.vec2

!
! 2) Test cross product of two vector fields
!
  resVec = vec1.cross.vec2

!
! 3) Tensor product of two vector fields
!
  T = vec1.tensor.vec2

  !...curl of a derived tensor field

  resVec = .curl.(3.14_dp*T+T)

  ! ...set field name for vtk output file name:
  resVec % field_name = 'Curl_volVectorField'

!
! 4) Magnitude squared of a tensor field
!
  phi = .magSq.T


  ! 6) Write output files in Paraview .vtu format
!+-----------------------------------------------------------------------------+

  ierr = write_volScalarField_field( phi )

  ierr = write_volVectorField_field( resVec )



  ! 7) Final touches.
!+-----------------------------------------------------------------------------+

   call say_goodbye

  end program

