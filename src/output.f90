module output
use types, only: dp
use parameters
use geometry,only: numCells,numInnerFaces,numBoundaryFaces,numFaces,numNodes,nomax,native_mesh_files,x,y,z,arx,ary,arz, &
                   read_line_faces_file_polyMesh
use tensor_fields
use utils, only: get_unit, i4_to_s_left

! 
! > Derived operators
!
interface operator(.visualize.)
   module procedure write_volScalarField_field
   module procedure write_volVectorField_field
end interface

public 

contains


 pure integer function paraview_cell_type(nnodes)
!
! Element type in Paraview corresponding to number of nodes.
!
 implicit none
 integer, intent(in) :: nnodes
 

  if(nnodes.eq.1) then 
   paraview_cell_type = 1
  elseif(nnodes.eq.2) then 
   paraview_cell_type =  3
  elseif(nnodes.eq.3) then 
   paraview_cell_type = 5
  elseif(nnodes.eq.4) then 
   paraview_cell_type = 9
  else
   paraview_cell_type = 0
  endif

 end function


 pure integer function paraview_ntype(NTYPE)
!
! Element type in Paraview corresponding to Cappuccino type of element given by NTYPE.
!
 implicit none
 integer, intent(in) :: NTYPE
 integer, parameter :: cappuccinoNTYPE1 = 3 !  Line
 integer, parameter :: cappuccinoNTYPE2 = 5 !  Tri
 integer, parameter :: cappuccinoNTYPE3 = 9 !  Quad
 integer, parameter :: cappuccinoNTYPE4 = 10 !  Tet
 integer, parameter :: cappuccinoNTYPE5 = 12 !  Hex
 integer, parameter :: cappuccinoNTYPE6 = 13 !  Prism
 integer, parameter :: cappuccinoNTYPE7 = 14 !  Pyramid
 

 if(NTYPE.eq.1) then ! -> Cappuccino line
   paraview_ntype = cappuccinoNTYPE1 
 elseif(NTYPE.eq.2) then ! -> Cappuccino Tri
   paraview_ntype = cappuccinoNTYPE2 
 elseif(NTYPE.eq.3) then ! -> Cappuccino Quad
   paraview_ntype = cappuccinoNTYPE3 
 elseif(NTYPE.eq.4) then ! -> Cappuccino Tet
   paraview_ntype = cappuccinoNTYPE4   
 elseif(NTYPE.eq.5) then ! -> Cappuccino Hex
   paraview_ntype = cappuccinoNTYPE5 
 elseif(NTYPE.eq.6) then ! -> Cappuccino Prism
   paraview_ntype = cappuccinoNTYPE6 
 elseif(NTYPE.eq.7) then ! -> Cappuccino Pyramid
   paraview_ntype = cappuccinoNTYPE7 
 else
   paraview_ntype = 0
 endif

 end function

 pure integer function noel(NTYPE)
!
! Number of nodes in element of NTYPE type of element
!
 implicit none
 integer, intent(in) :: NTYPE
 integer, parameter :: noel_NTYPE1 = 2 ! no. of nodes in element-here Line
 integer, parameter :: noel_NTYPE2 = 3 ! no. of nodes in element-here Tri
 integer, parameter :: noel_NTYPE3 = 4 ! no. of nodes in element-here Quad
 integer, parameter :: noel_NTYPE4 = 4 ! no. of nodes in element-here Tet
 integer, parameter :: noel_NTYPE5 = 8 ! no. of nodes in element-here Hex
 integer, parameter :: noel_NTYPE6 = 6 ! no. of nodes in element-here Prism
 integer, parameter :: noel_NTYPE7 = 5 ! no. of nodes in element-here Pyramid
 

 if(NTYPE.eq.1) then ! -> Cappuccino line
   noel = noel_NTYPE1 
 elseif(NTYPE.eq.2) then ! -> Cappuccino Tri
   noel = noel_NTYPE2 
 elseif(NTYPE.eq.3) then ! -> Cappuccino Quad
   noel = noel_NTYPE3 
 elseif(NTYPE.eq.4) then ! -> Cappuccino Tet
   noel = noel_NTYPE4 
 elseif(NTYPE.eq.5) then ! -> Cappuccino Hex
   noel = noel_NTYPE5 
 elseif(NTYPE.eq.6) then ! -> Cappuccino Prism
   noel = noel_NTYPE6 
 elseif(NTYPE.eq.7) then ! -> Cappuccino Pyramid
   noel = noel_NTYPE7 
 else
   noel = 0
 endif

 end function



function write_volScalarField_field(phi) result(ierr)
!
! Wrapper around subroutine vtu_write for VolScalarField
!
  implicit none
  type(volScalarField), intent(in) :: phi
  integer :: ierr

  integer :: output_unit

  call get_unit( output_unit )

  open(unit = output_unit, file = trim(phi%field_name)//'_field_data.vtu', status = 'replace', iostat = ierr)

  call vtu_write ( output_unit, trim(phi%field_name)//'_field_data', phi%mag )

  close (  unit = output_unit )

end function


function write_volVectorField_field(v) result(ierr)
!
! Wrapper around subroutine vtu_write for VolScalarField
!
  implicit none
  type(volVectorField), intent(in) :: v
  integer :: ierr

  integer :: output_unit

  call get_unit( output_unit )

  open(unit = output_unit, file = trim(v%field_name)//'_field_data.vtu', status = 'replace', iostat = ierr)

  call vtu_write_vector_field ( output_unit, trim(v%field_name)//'_field_data', v%x, v%y, v%z )

  close (  unit = output_unit )

end function



subroutine vtu_write ( output_unit, scalar_name, scalar_field )
!
! Writes scalar field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit
  character ( len = * ), intent(in) :: scalar_name
  real(dp), dimension(numCells), intent(in) :: scalar_field

  character ( len = 20 ) node_num_string
  character ( len = 20 ) cells_num_string
  ! character ( len = 1 ) ch

  integer :: i,k
  integer :: icell
  ! integer :: iface
  integer :: ntype
  integer :: offset
  integer :: cells_file
  ! integer :: faces_file
  ! integer :: nnodes
  integer, dimension(8) :: node
  ! integer, dimension(numBoundaryFaces) :: num_nodes

!
! > Header
!

  call i4_to_s_left ( numNodes, node_num_string )
  call i4_to_s_left ( numCells, cells_num_string )
  ! call i4_to_s_left ( numCells+numBoundaryFaces, cells_num_string )

  call get_unit( cells_file )
  open( unit = cells_file, file='polyMesh/cells' )
  rewind cells_file

  ! call get_unit( faces_file )
  ! open( unit = faces_file, file='polyMesh/faces' )
  ! rewind faces_file
  ! ! Fast forward to place where boundary faces start
  !   ! First polyMesh header
  !   do i=1,18
  !     read(faces_file,*) ch
  !   end do
  !   ! Then over Inner faces
  !   do i=1,numInnerFaces
  !     read(faces_file,*) ch
  !   end do

  write ( output_unit, '(a)' )    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
  write ( output_unit, '(2x,a)' ) '<UnstructuredGrid>'
  write ( output_unit, '(4x,5a)' ) &
  '<Piece NumberOfPoints="',trim( node_num_string ),'" NumberOfCells="',trim( cells_num_string ),'">'


! ! <Scalars in nodes>
!  write ( output_unit, '(6x,a)' ) '<PointData Scalars="scalars">'
!  write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" Name="',scalar_name,'" Format="ascii">'

! !   <Scalar field data>
!    do i=1,numNodes
!      write( output_unit, '(es11.4)') scalar_field(i) 
!    enddo
! !   </Scalar field data>

!  write ( output_unit, '(8x,a)' ) '</DataArray>'
!  write ( output_unit, '(6x,a)' ) '</PointData>'
! ! </Scalars in nodes>

! <Scalars in cell-centers>
  write ( output_unit, '(6x,a)' ) '<CellData Scalars="scalars">'
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',scalar_name,'" Format="ascii">'

!
! > Scalars in cell-centers > write scalar data
!
    do icell=1,numCells
      write( output_unit, '(10x,es11.4)') scalar_field(icell) 
    enddo
    ! do iface = 1,numBoundaryFaces
    !   i = numCells+iface
    !   write( output_unit, '(10x,es11.4)') scalar_field(i)
    ! enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</CellData>'
! </Scalars in cell-centers>

!
! > Mesh data
!
  write ( output_unit, '(6x,a)' ) '<Points>'
  write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

!
! > Mesh data > Nodal coordinates
!
  do i=1,numNodes
    write ( output_unit, '(10x,3es11.4)' ) x(i),y(i),z(i)
  enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</Points>'
!
! > Mesh data > Connectivity
!
  write ( output_unit, '(6x,a)' ) '<Cells>'

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
    do icell=1,numCells
      read( cells_file, * ) ntype,(node(k), k=1,noel(ntype))
      ! Note, Paraview starts counting from 0, we in Fortran start with 1, therefore: node(k)-1
      write( output_unit, '(10x,4i8:(4i8:))') (node(k)-1, k=1,noel(ntype)) 
    enddo
    ! do iface=1,numBoundaryFaces
    !   ! Read line in 'faces' file
    !   if (native_mesh_files) then
    !     read( faces_file, * ) nnodes, (node(k), k=1,nnodes)
    !   else ! OpenFOAM polyMesh
    !     call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)
    !   endif
    !   ! Note, Paraview starts counting from 0, so does native polyMesh, but in 
    !   ! Cappuccino we have to substract 1 as above for cells.
    !   if (native_mesh_files) then
    !     write( output_unit, '(10x,4i8:(4i8:))') ( node(k)-1, k=1,nnodes ) 
    !   else
    !     write( output_unit, '(10x,4i8:(4i8:))') ( node(k), k=1,nnodes ) 
    !   endif
    !   num_nodes(iface) = nnodes
    ! enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  rewind cells_file

  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'
    do icell=1,numCells
      read( cells_file, '(i1)' ) ntype
      offset = offset+noel(ntype)
      write( output_unit, '(i12)') offset 
    enddo
    ! do iface=1,numBoundaryFaces
    !   offset = offset+num_nodes(iface)
    !   write( output_unit, '(i12)') offset 
    ! enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

  rewind cells_file
!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'
    do icell=1,numCells
      read( cells_file, '(i1)' ) ntype
      write( output_unit, '(i12)') paraview_ntype(ntype) 
    enddo
    ! do iface=1,numBoundaryFaces
    !   write( output_unit, '(i12)') paraview_cell_type( num_nodes(iface) ) 
    ! enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(6x,a)' ) '</Cells>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
  write ( output_unit, '(a)' ) '</VTKFile>'

  close ( cells_file )

end subroutine vtu_write




subroutine vtu_write_vector_field ( output_unit, field_name, u, v, w )
!
! Writes scalar field data to Paraview XML, unstructured, ".vtu" file.
!

  implicit none
  integer, intent(in) :: output_unit
  character ( len = * ), intent(in) :: field_name
  real(dp), dimension(numCells), intent(in) :: u, v, w

  character ( len = 20 ) node_num_string
  character ( len = 20 ) cells_num_string
  ! character ( len = 1 ) ch

  integer :: i,k
  integer :: icell
  ! integer :: iface
  integer :: ntype
  integer :: offset
  integer :: cells_file
  ! integer :: faces_file
  ! integer :: nnodes
  integer, dimension(8) :: node
  ! integer, dimension(numBoundaryFaces) :: num_nodes

!
! > Header
!

  call i4_to_s_left ( numNodes, node_num_string )
  call i4_to_s_left ( numCells, cells_num_string )
  ! call i4_to_s_left ( numCells+numBoundaryFaces, cells_num_string )

  call get_unit( cells_file )
  open( unit = cells_file, file='polyMesh/cells' )
  rewind cells_file

  ! call get_unit( faces_file )
  ! open( unit = faces_file, file='polyMesh/faces' )
  ! rewind faces_file
  ! ! Fast forward to place where boundary faces start
  !   ! First polyMesh header
  !   do i=1,18
  !     read(faces_file,*) ch
  !   end do
  !   ! Then over Inner faces
  !   do i=1,numInnerFaces
  !     read(faces_file,*) ch
  !   end do

  write ( output_unit, '(a)' )    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
  write ( output_unit, '(2x,a)' ) '<UnstructuredGrid>'
  write ( output_unit, '(4x,5a)' ) &
  '<Piece NumberOfPoints="',trim( node_num_string ),'" NumberOfCells="',trim( cells_num_string ),'">'


! <Scalars in cell-centers>
  write ( output_unit, '(6x,a)' ) '<CellData Scalars="scalars">'
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',trim( field_name ),'" NumberOfComponents="3" Format="ascii">'

!
! > Scalars in cell-centers > write scalar data
!
    do icell=1,numCells
      write( output_unit, '(10x,3(1x,es11.4))') u(icell), v(icell), w(icell)
    enddo
    ! do iface = 1,numBoundaryFaces
    !   i = numCells+iface
    !   write( output_unit, '(10x,3(1x,es11.4))') u(i), v(i), w(i)
    ! enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</CellData>'
! </Scalars in cell-centers>

!
! > Mesh data
!
  write ( output_unit, '(6x,a)' ) '<Points>'
  write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

!
! > Mesh data > Nodal coordinates
!
  do i=1,numNodes
    write ( output_unit, '(10x,3es11.4)' ) x(i),y(i),z(i)
  enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</Points>'
!
! > Mesh data > Connectivity
!
  write ( output_unit, '(6x,a)' ) '<Cells>'

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
    do icell=1,numCells
      read( cells_file, * ) ntype,(node(k), k=1,noel(ntype))
      ! Note, Paraview starts counting from 0, we in Fortran start with 1, therefore: node(k)-1
      write( output_unit, '(10x,4i8:(4i8:))') (node(k)-1, k=1,noel(ntype)) 
    enddo
    ! do iface=1,numBoundaryFaces
    !   ! Read line in 'faces' file
    !   if (native_mesh_files) then
    !     read( faces_file, * ) nnodes, (node(k), k=1,nnodes)
    !   else ! OpenFOAM polyMesh
    !     call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)
    !   endif
    !   ! Note, Paraview starts counting from 0, so does native polyMesh, but in 
    !   ! Cappuccino we have to substract 1 as above for cells.
    !   if (native_mesh_files) then
    !     write( output_unit, '(10x,4i8:(4i8:))') ( node(k)-1, k=1,nnodes ) 
    !   else
    !     write( output_unit, '(10x,4i8:(4i8:))') (node(k), k=1,nnodes ) 
    !   endif
    !   num_nodes(iface) = nnodes
    ! enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  rewind cells_file
  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'
    do icell=1,numCells
      read( cells_file, '(i1)' ) ntype
      offset = offset+noel(ntype)
      write( output_unit, '(i12)') offset 
    enddo
    ! do iface=1,numBoundaryFaces
    !   offset = offset+num_nodes(iface)
    !   write( output_unit, '(i12)') offset 
    ! enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

  rewind cells_file
!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'
    do icell=1,numCells
      read( cells_file, '(i1)' ) ntype
      write( output_unit, '(i12)') paraview_ntype(ntype) 
    enddo
    ! do iface=1,numBoundaryFaces
    !   write( output_unit, '(i12)') paraview_cell_type( num_nodes(iface) ) 
    ! enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(6x,a)' ) '</Cells>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
  write ( output_unit, '(a)' ) '</VTKFile>'

  close ( cells_file )

end subroutine vtu_write_vector_field





subroutine vtu_write_mesh ( output_unit )
!
! Writes scalar field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit

  character ( len = 20 ) node_num_string
  character ( len = 20 ) cells_num_string

  integer :: i,k
  integer :: icell
  integer :: ntype
  integer :: offset
  integer :: cells_file

  integer, dimension(8) :: node

!
! > Header
!

  call i4_to_s_left ( numNodes, node_num_string )
  call i4_to_s_left ( numCells, cells_num_string )

  call get_unit( cells_file )
  open( unit = cells_file, file='polyMesh/cells' )
  rewind cells_file

  write ( output_unit, '(a)' )    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
  write ( output_unit, '(2x,a)' ) '<UnstructuredGrid>'
  write ( output_unit, '(4x,5a)' ) &
  '<Piece NumberOfPoints="',trim( node_num_string ),'" NumberOfCells="',trim( cells_num_string ),'">'

!
! > Mesh data
!
  write ( output_unit, '(6x,a)' ) '<Points>'
  write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

!
! > Mesh data > Nodal coordinates
!
  do i=1,numNodes
    write ( output_unit, '(10x,3es11.4)' ) x(i),y(i),z(i)
  enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</Points>'
!
! > Mesh data > Connectivity
!
  write ( output_unit, '(6x,a)' ) '<Cells>'

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
    do icell=1,numCells
      read( cells_file, '(i1,1x,4i8:(4i8:))' ) ntype,(node(k), k=1,noel(ntype))
      ! Note, Paraview starts counting from 0, we in Fortran start with 1, therefore: node(k)-1
      write( output_unit, '(4i8:(4i8:))') (node(k)-1, k=1,noel(ntype)) 
    enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  rewind cells_file
  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'
    do icell=1,numCells
      read( cells_file, '(i1)' ) ntype
      offset = offset+noel(ntype)
      write( output_unit, '(i12)') offset 
    enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

  rewind cells_file
!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'
    do icell=1,numCells
      read( cells_file, '(i1)' ) ntype
      write( output_unit, '(i12)') paraview_ntype(ntype) 
    enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(6x,a)' ) '</Cells>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
  write ( output_unit, '(a)' ) '</VTKFile>'

  close ( cells_file )

end subroutine vtu_write_mesh



subroutine vtp_write ( output_unit, scalar_name, scalar_field )
!
! Writes scalar field data to Paraview POLYDATA format, ".vtp" file.
!
  implicit none

  integer, intent(in) :: output_unit
  character ( len = * ), intent(in) :: scalar_name
  real(dp), dimension(numCells), intent(in) :: scalar_field

  character ( len = 20 ) node_num_string
  character ( len = 20 ) cells_num_string
  character ( len = 20)  faces_num_string
  character ( len = 1) ch

  integer :: i,k
  integer :: npts
  integer :: offset
  integer :: faces_file
  integer :: nnodes 
  integer, dimension(nomax) :: node

  real(dp) :: are, nx, ny, nz

!
! > Header
!

  call i4_to_s_left ( numNodes, node_num_string )
  call i4_to_s_left ( numCells, cells_num_string )
  call i4_to_s_left ( numFaces, faces_num_string )

  call get_unit( faces_file )
  open( unit = faces_file, file='polyMesh/faces' )
  rewind faces_file

  do i=1,100
    read(faces_file,'(a)') ch
    if (ch == "(") exit
  end do

  write ( output_unit, '(a)' )    '<?xml version="1.0"?>'
  write ( output_unit, '(a)' )    '<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">'
  write ( output_unit, '(2x,a)' ) '<PolyData>'
  write ( output_unit, '(4x,5a)' ) '<Piece NumberOfPoints="',trim( node_num_string ), &
  '" NumberOfVerts="0" NumberOfLines="0" NumberOfStrips="0" NumberOfPolys="',trim( faces_num_string ),'">'


!
! > Mesh data
!
  write ( output_unit, '(6x,a)' ) '<Points>'
  write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'

!
! > Mesh data > Nodal coordinates
!
  do i=1,numNodes
    write ( output_unit, '(10x,3es11.4)' ) x(i),y(i),z(i)
  enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</Points>'


! <Scalars in nodes>
 write ( output_unit, '(6x,a)' ) '<PointData Scalars="my_scalars">'
 write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" Name="',scalar_name,'" Format="ascii">'

!   <Scalar field data>
    do i=1,numNodes
     write( output_unit, '(es11.4)') scalar_field(i) 
    enddo
!   </Scalar field data>

 write ( output_unit, '(8x,a)' ) '</DataArray>'
 write ( output_unit, '(6x,a)' ) '</PointData>'
! </Scalars in nodes>


! <Scalars in cell-centers>
  write ( output_unit, '(6x,a)' ) '<CellData Scalars="face_scalars" Normals="cell_normals">'
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',scalar_name,'" Format="ascii">'

!
! > Scalars in cell-centers > write scalar data
!
    do i=1,numFaces
      write( output_unit, '(10x,es11.4)') scalar_field(i) 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(8x,a)' ) '<DataArray type="Float32" Name="cell_normals" NumberOfComponents="3" Format="ascii">'

!
! > Scalars in cell-centers > write cell normals
!
    do i=1,numFaces

      are = sqrt(arx(i)**2+ary(i)**2+arz(i)**2)

      nx = arx(i)/are
      ny = ary(i)/are
      nz = arz(i)/are

      write( output_unit, '(10x,3es11.4)') nx,ny,nz 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'
  write ( output_unit, '(6x,a)' ) '</CellData>'
! </Scalars in cell-centers>


!
! > Mesh data > Connectivity
!
  write ( output_unit, '(6x,a)' ) '<Polys>'

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="connectivity" Format="ascii">'


  do i=1,numFaces

    node(:) = 0

    ! Read line in 'faces' file
    if (native_mesh_files) then
      read( faces_file, * ) nnodes, (node(k), k=1,nnodes)
    else ! OpenFOAM polyMesh
      call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)
    endif

    ! Note, Paraview starts counting from 0, we in Fortran start with 1, therefore: node(k)-1
    write( output_unit, '(10x,12i8:(12i8:))') (node(k)-1, k=1,nnodes) 

  enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  rewind faces_file
  do i=1,100
    read(faces_file,'(a)') ch
    if (ch == "(") exit
  end do

  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'
    do i=1,numFaces
      read( faces_file, '(i1)' ) npts
      offset = offset + npts
      write( output_unit, '(i12)') offset 
    enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'


  write ( output_unit, '(6x,a)' ) '</Polys>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</PolyData>'
  write ( output_unit, '(a)' ) '</VTKFile>'

  close ( faces_file )

end subroutine vtp_write


end module
