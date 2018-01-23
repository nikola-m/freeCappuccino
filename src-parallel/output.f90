module output
use types, only: dp
use parameters
use geometry
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

 if(NTYPE.eq.1) then ! -> Cappuccino line
   paraview_ntype = 3 
 elseif(NTYPE.eq.2) then ! -> Cappuccino Tri
   paraview_ntype = 5 
 elseif(NTYPE.eq.3) then ! -> Cappuccino Quad
   paraview_ntype = 9 
 elseif(NTYPE.eq.4) then ! -> Cappuccino Tet
   paraview_ntype = 10   
 elseif(NTYPE.eq.5) then ! -> Cappuccino Hex
   paraview_ntype = 12 
 elseif(NTYPE.eq.6) then ! -> Cappuccino Prism
   paraview_ntype = 13 
 elseif(NTYPE.eq.7) then ! -> Cappuccino Pyramid
   paraview_ntype = 14 
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


 if(NTYPE.eq.12) then ! -> Paraview Hex
   noel = 8  
 elseif(NTYPE.eq.13) then ! -> Paraview Prism
   noel = 6 
 elseif(NTYPE.eq.14) then ! -> Paraview Pyramid
   noel = 5
 elseif(NTYPE.eq.10) then ! -> Paraview Tet
   noel = 4
 elseif(NTYPE.eq.9) then ! -> Paraview Quad
   noel = 4
 elseif(NTYPE.eq.5) then ! -> Paraview Tri
   noel = 3
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

  call vtu_write_scalar_field ( output_unit, trim(phi%field_name)//'_field_data', phi%mag )

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



subroutine vtu_write_XML_header ( output_unit )
!
! Writes header data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit

  integer :: icell

  character ( len = 20 ) node_num_string
  character ( len = 20 ) cells_num_string

!
! > Header
!

  call i4_to_s_left ( numNodes, node_num_string )
  call i4_to_s_left ( numCells+numBoundaryFaces, cells_num_string )


  write ( output_unit, '(a)' )    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
  write ( output_unit, '(2x,a)' ) '<UnstructuredGrid>'
  write ( output_unit, '(4x,5a)' ) &
  '<Piece NumberOfPoints="',trim( node_num_string ),'" NumberOfCells="',trim( cells_num_string ),'">'

! > Open XML tag for recording data in cell centers
  write ( output_unit, '(6x,a)' ) '<CellData>'

  call i4_to_s_left ( numCells+numBoundaryFaces-1, cells_num_string )

  write ( output_unit, '(8x,a)' ) &
  & '<DataArray type="Int32" Name="cellID" format="ascii" RangeMin="0" RangeMax="'//trim( cells_num_string )//'">'

    do icell=1,numCells+numBoundaryFaces
      write( output_unit, '(10x,i8)') icell
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

end subroutine


subroutine vtu_write_XML_meshdata ( output_unit )
!
! Writes unstructured mesh data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit

  integer :: i,k
  integer :: icell
  integer :: ntype
  integer :: offset
  integer :: cells_file
  integer, dimension(8) :: node
  character( len = 5) :: nproc_char

  ! nproc_char <- myid zapisan levo u vidu stringa.
  call i4_to_s_left ( myid, nproc_char )

  call get_unit( cells_file )
  open( unit = cells_file, file='processor'//trim(nproc_char)//'/constant/polyMesh/cells' )
  rewind cells_file


! > End of recording data in cell centers
  write ( output_unit, '(6x,a)' ) '</CellData>'

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

    do icell=1,numCells+numBoundaryFaces
      read( cells_file, * ) ntype,(node(k), k=1,noel(ntype))
      write( output_unit, '(10x,4i8:(4i8:))') (node(k), k=1,noel(ntype)) 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  rewind cells_file

  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'

    do icell=1,numCells+numBoundaryFaces
      read( cells_file, * ) ntype
      offset = offset+noel(ntype)
      write( output_unit, '(10x,i12)') offset 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  rewind cells_file
!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'

    do icell=1,numCells+numBoundaryFaces
      read( cells_file, * ) ntype
      write( output_unit, '(10x,i2)') ntype 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(6x,a)' ) '</Cells>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
  write ( output_unit, '(a)' ) '</VTKFile>'

  close ( cells_file )

end subroutine



subroutine vtu_write_XML_scalar_field ( output_unit, scalar_name, scalar_field )
!
! Writes scalar field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit
  character ( len = * ), intent(in) :: scalar_name
  real(dp), dimension(numCells), intent(in) :: scalar_field

  integer :: icell
  integer :: i,k
  integer :: nfaces
  integer :: startFace
  integer :: boundary_file
  integer :: numBoundaries
  integer :: numinl 
  integer :: numout
  integer :: numsym
  integer :: numwal
  integer :: numpru
  integer :: numwali
  integer :: numwala
  integer :: numwalf
  character(len=10) :: bctype
  character(len=80) :: line_string
  character( len = 5) :: nproc_char


  ! nproc_char <- myid zapisan levo u vidu stringa.
  call i4_to_s_left ( myid, nproc_char )

! <Scalars in cell-centers>
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',scalar_name,'" Format="ascii">'

!
! > Scalars in cell-centers > write scalar data
!
    do icell=1,numCells
      write( output_unit, '(10x,es11.4)') scalar_field(icell) 
    enddo

!
! > Scalars in boundary faces > write scalar data
!

  numinl = 0 
  numout = 0
  numsym = 0
  numwal = 0
  numpru = 0
  numwali = 0
  numwala = 0
  numwalf = 0

  call get_unit( boundary_file )
  open( unit = boundary_file, file='processor'//trim(nproc_char)//'/constant/polyMesh/boundary' )
  rewind boundary_file

  ! Number of rows in the file excluding #comment in header
  call file_row_count ( boundary_file, numBoundaries )

  read(boundary_file,'(a)') line_string

  do k=1,numBoundaries
  read(boundary_file,*) bctype,nfaces,startFace

    select case ( adjustl(bctype) )

      case ('inlet')

        do i=1,nfaces
          icell = iInletStart + numinl + i
          write( output_unit, '(10x,es11.4)') scalar_field(icell)
        enddo

        numinl = numinl + nfaces

      case ('outlet')

        do i=1,nfaces
          icell = iOutletStart + numout + i
          write( output_unit, '(10x,es11.4)') scalar_field(icell)
        enddo

        numout = numout + nfaces

      case ('symmetry')

        do i=1,nfaces
          icell = iSymmetryStart + numsym + i
          write( output_unit, '(10x,es11.4)') scalar_field(icell)
        enddo

        numsym = numsym + nfaces

      case ('wall')

        do i=1,nfaces
          icell = iWallStart + numwal + i
          print*, icell, iWallStart, numwal, i
          write( output_unit, '(10x,es11.4)') scalar_field(icell)
        enddo

        numwal = numwal + nfaces

      case ('wallIsoth')

        do i=1,nfaces
          icell = iWallStart + numwal + i
          write( output_unit, '(10x,es11.4)') scalar_field(icell)
        enddo

        numwal = numwal + nfaces

      case ('wallAdiab')

        do i=1,nfaces
          icell = iWallStart + numwal + i
          write( output_unit, '(10x,es11.4)') scalar_field(icell)
        enddo

        numwal = numwal + nfaces

      case ('wallQFlux')

        do i=1,nfaces
          icell = iWallStart + numwal + i
          write( output_unit, '(10x,es11.4)') scalar_field(icell)
        enddo

        numwal = numwal + nfaces

      case ('prOutlet')

        do i=1,nfaces
          icell = iPressOutletStart + numpru + i
          write( output_unit, '(10x,es11.4)') scalar_field(icell)
        enddo

        numpru = numpru + nfaces

      case default
        write(*,*) "Output.f90: Non-existing boundary type in polymesh/boundary file!"
        stop

    end select

  enddo

  close ( boundary_file )

  write ( output_unit, '(8x,a)' ) '</DataArray>'
! </Scalars in cell-centers>

end subroutine vtu_write_XML_scalar_field



subroutine vtu_write_XML_vector_field ( output_unit, field_name, u, v, w )
!
! Writes vector field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit
  character ( len = * ), intent(in) :: field_name
  real(dp), dimension(numCells), intent(in) :: u, v, w

  integer :: icell
  integer :: i,k
  integer :: nfaces
  integer :: startFace
  integer :: boundary_file
  integer :: numBoundaries
  integer :: numinl 
  integer :: numout
  integer :: numsym
  integer :: numwal
  integer :: numpru
  integer :: numwali
  integer :: numwala
  integer :: numwalf
  character(len=10) :: bctype
  character(len=80) :: line_string
  character( len = 5) :: nproc_char

  ! nproc_char <- myid zapisan levo u vidu stringa.
  call i4_to_s_left ( myid, nproc_char )

! <Scalars in cell-centers>
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',trim( field_name ),'" NumberOfComponents="3" Format="ascii">'

!
! > Scalars in cell-centers > write scalar data
!
    do icell=1,numCells
      write( output_unit, '(10x,3(1x,es11.4))') u(icell), v(icell), w(icell)
    enddo

!
! > Scalars in boundary faces > write scalar data
!

  numinl = 0 
  numout = 0
  numsym = 0
  numwal = 0
  numpru = 0
  numwali = 0
  numwala = 0
  numwalf = 0

  call get_unit( boundary_file )
  open( unit = boundary_file, file='processor'//trim(nproc_char)//'/constant/polyMesh/boundary' )
  rewind boundary_file

  ! Number of rows in the file excluding #comment in header
  call file_row_count ( boundary_file, numBoundaries )

  read(boundary_file,'(a)') line_string

  do k=1,numBoundaries
  read(boundary_file,*) bctype,nfaces,startFace

    select case ( adjustl(bctype) )

      case ('inlet')

        do i=1,nfaces
          icell = iInletStart + numinl + i
          write( output_unit, '(10x,es11.4)') u(icell), v(icell), w(icell)
        enddo

        numinl = numinl + nfaces

      case ('outlet')

        do i=1,nfaces
          icell = iOutletStart + numout + i
          write( output_unit, '(10x,es11.4)') u(icell), v(icell), w(icell)
        enddo

        numout = numout + nfaces

      case ('symmetry')

        do i=1,nfaces
          icell = iSymmetryStart + numsym + i
          write( output_unit, '(10x,es11.4)') u(icell), v(icell), w(icell)
        enddo

        numsym = numsym + nfaces

      case ('wall')

        do i=1,nfaces
          icell = iWallStart + numwal + i
          write( output_unit, '(10x,es11.4)') u(icell), v(icell), w(icell)
        enddo

        numwal = numwal + nfaces

      case ('wallIsoth')

        do i=1,nfaces
          icell = iWallStart + numwal + i
          write( output_unit, '(10x,es11.4)') u(icell), v(icell), w(icell)
        enddo

        numwal = numwal + nfaces

      case ('wallAdiab')

        do i=1,nfaces
          icell = iWallStart + numwal + i
          write( output_unit, '(10x,es11.4)') u(icell), v(icell), w(icell)
        enddo

        numwal = numwal + nfaces

      case ('wallQFlux')

        do i=1,nfaces
          icell = iWallStart + numwal + i
          write( output_unit, '(10x,es11.4)') u(icell), v(icell), w(icell)
        enddo

        numwal = numwal + nfaces

      case ('prOutlet')

        do i=1,nfaces
          icell = iPressOutletStart + numpru + i
          write( output_unit, '(10x,es11.4)') u(icell), v(icell), w(icell)
        enddo

        numpru = numpru + nfaces
        
      case default
        write(*,*) "Output.f90: Non-existing boundary type in polymesh/boundary file!"
        stop

    end select

  enddo

  close ( boundary_file )


  write ( output_unit, '(8x,a)' ) '</DataArray>'
! </Scalars in cell-centers>

end subroutine vtu_write_XML_vector_field



subroutine vtu_write_scalar_field ( output_unit, scalar_name, scalar_field )
!
! Writes scalar field data to Paraview XML, unstructured, ".vtu" file.
!
  implicit none

  integer, intent(in) :: output_unit
  character ( len = * ), intent(in) :: scalar_name
  real(dp), dimension(numCells), intent(in) :: scalar_field

  character ( len = 20 ) node_num_string
  character ( len = 20 ) cells_num_string
  character( len = 5) :: nproc_char

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
  call i4_to_s_left ( myid, nproc_char )

  call get_unit( cells_file )
  open( unit = cells_file, file='processor'//trim(nproc_char)//'/constant/polyMesh/cells' )
  rewind cells_file


  write ( output_unit, '(a)' )    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">'
  write ( output_unit, '(2x,a)' ) '<UnstructuredGrid>'
  write ( output_unit, '(4x,5a)' ) &
  '<Piece NumberOfPoints="',trim( node_num_string ),'" NumberOfCells="',trim( cells_num_string ),'">'


! <Scalars in cell-centers>
  write ( output_unit, '(6x,a)' ) '<CellData Scalars="scalars">'
  write ( output_unit, '(8x,3a)' ) '<DataArray type="Float32" Name="',scalar_name,'" Format="ascii">'

!
! > Scalars in cell-centers > write scalar data
!
    do icell=1,numCells
      write( output_unit, '(10x,es11.4)') scalar_field(icell) 
    enddo

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
      write( output_unit, '(10x,4i8:(4i8:))') (node(k), k=1,noel(ntype)) 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  rewind cells_file

  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype
      offset = offset+noel(ntype)
      write( output_unit, '(10x,i12)') offset 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  rewind cells_file
!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype
      write( output_unit, '(10x,i2)') ntype 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(6x,a)' ) '</Cells>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
  write ( output_unit, '(a)' ) '</VTKFile>'

  close ( cells_file )

end subroutine vtu_write_scalar_field




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
  character( len = 5) :: nproc_char

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
  call i4_to_s_left ( myid, nproc_char )

  call get_unit( cells_file )
  open( unit = cells_file, file='processor'//trim(nproc_char)//'/constant/polyMesh/cells' )
  rewind cells_file

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
      write( output_unit, '(10x,4i8:(4i8:))') (node(k), k=1,noel(ntype)) 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  rewind cells_file
  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype
      offset = offset+noel(ntype)
      write( output_unit, '(10x,i12)') offset 
    enddo

  write ( output_unit, '(8x,a)' ) '</DataArray>'

  rewind cells_file
!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'

    do icell=1,numCells
      read( cells_file, * ) ntype
      write( output_unit, '(10x,i2)') ntype
    enddo

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
  character( len = 5) :: nproc_char

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
  call i4_to_s_left ( myid, nproc_char )

  call get_unit( cells_file )
  open( unit = cells_file, file='processor'//trim(nproc_char)//'/constant/polyMesh/cells' )
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
      write( output_unit, '(4i8:(4i8:))') (node(k), k=1,noel(ntype)) 
    enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

!
! > Mesh data > Offsets
!
  rewind cells_file
  offset = 0

  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="offsets" Format="ascii">'
    do icell=1,numCells
      read( cells_file, * ) ntype
      offset = offset+noel(ntype)
      write( output_unit, '(10x,i12)') offset 
    enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

  rewind cells_file
!
! > Mesh data > Types
!
  write ( output_unit, '(8x,a)' ) '<DataArray type="Int32" Name="types" Format="ascii">'
    do icell=1,numCells
      read( cells_file, * ) ntype
      write( output_unit, '(10x,i2)') ntype 
    enddo
  write ( output_unit, '(8x,a)' ) '</DataArray>'

  write ( output_unit, '(6x,a)' ) '</Cells>'
  write ( output_unit, '(4x,a)' ) '</Piece>'
  write ( output_unit, '(2x,a)' ) '</UnstructuredGrid>'
  write ( output_unit, '(a)' ) '</VTKFile>'

  close ( cells_file )

end subroutine vtu_write_mesh


end module