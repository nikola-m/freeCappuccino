!*******************************************************************************
program writeCells
!
!  Description:
!    Read mesh in polyMesh format and creates messy 'cells' file, which is 
!    to be used in output routines to write postprocessing files 
!    for Paraview in .vtu unstructured mesh format
!    The point is to get the file with entries:
!      ...
!      NNODES NODE_INDX_1 NODE_INDX_2 ... NODE_INDX_NNODES
!      ...
!    Where node indices may repeat, but the number of nodes in the line
!    must be NNODES.
!    Input mesh is described in files similar to polyMesh format.
!    Mesh files are 'points', 'faces', 'owner, 'neighbour', 'boundary'
!
!   
!
!  Date:
!    13/10/2017
!
!  Author:
!    Nikola Mirkov nmirkov@vinca.rs
!
!*******************************************************************************
!
  use cells_module
  use qsort_c_module

  implicit none


  ! Logical defining are we reading native Cappuccino file format or OpenFOAM polymesh format.
  logical :: native_mesh_files

  ! Mesh file units
  integer :: cells_file, faces_file, owner_file, neighbour_file

  character(len=1) :: ch

  integer, parameter :: nomax = 8 ! Change this if necessary, for size.
  integer :: iface
  integer :: i,k
  integer :: inp,inn
  integer :: nnodes
  integer :: numFaces, numInnerFaces, numBoundaryFaces
  integer :: node(nomax) 
  integer, dimension(:), allocatable :: cInd

!
! > OPEN polyMesh format files: 'points', 'faces', 'owner', 'neighbour'.
!

  call get_unit( faces_file )
  open( unit = faces_file, file='polyMesh/faces' )
  rewind faces_file

  call get_unit( owner_file )
  open( unit = owner_file, file='polyMesh/owner' )
  rewind owner_file

  call get_unit( neighbour_file )
  open( unit = neighbour_file, file='polyMesh/neighbour' )
  rewind neighbour_file

  call get_unit( cells_file )
  open( unit = cells_file, file='polyMesh/cells' )
  rewind cells_file

  native_mesh_files = .true.

  ! First characters in OpenFOAM polyMesh files is comment '/*'
  read(owner_file,'(a)') ch
  backspace(owner_file)
  if(ch == '/')  native_mesh_files = .false.


  !
  ! Fast forward files to place where entries are
  !

  ! 'owner' file
  do i=1,17
      read(owner_file,*) ch
  end do
  read(owner_file,*) numFaces
  read(owner_file,*) ch ! reads "("


  ! 'neighbour' file
  do i=1,17
    read(neighbour_file,*) ch
  end do
  read(neighbour_file,*) numInnerFaces
  read(neighbour_file,*) ch ! reads "("


  numBoundaryFaces = numFaces - numInnerFaces


  allocate ( cInd(2*numInnerFaces+numBoundaryFaces) )


  ! 'faces' file
  do i=1,18
    read(faces_file,*) ch
  end do

  i = 0

  do iface=1,numFaces

    read(owner_file,*) inp
    read(owner_file,*) inn

    node = 0

    ! Read line in 'faces' file
    if (native_mesh_files) then
      read( faces_file, * ) nnodes, (node(k), k=1,nnodes)
    else ! OpenFOAM polyMesh
      call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)
    endif

    write(cells_file,*) inp, nnodes, node(1:nnodes)
    write(cells_file,*) inn, nnodes, node(1:nnodes)

    i = i+1
    cInd(i) = inp
    
    i = i+1
    cInd(i) = inn


  enddo

 close(cells_file)
 close(neighbour_file)
 close(owner_file)
 close(faces_file)

end program