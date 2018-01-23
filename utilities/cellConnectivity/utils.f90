module utils

 implicit none

 
 public


 contains
 

subroutine read_line_faces_file_polyMesh(faces_file,nn,nod,nmax)
  implicit none
  integer, intent(in) :: faces_file
  integer, intent(in) :: nmax
  integer, intent(out) :: nn
  integer, dimension(nmax), intent(out) :: nod
  integer :: j,m,n
  character(len=15) :: char_string,char_string2

    nn = 0
    nod(:) = 0

    ! Read how many nodes in face
    read(faces_file,'(a)') char_string
    read(char_string(1:1),*)j
    backspace(faces_file)

      read(faces_file,*) char_string,nod(2:j-1),char_string2

      ! Char to double conversion:
      read(char_string(1:1),*)j                ! number of vertrices
      read(char_string(3:),*) m                ! first vertex
      read(char_string2(1:len_trim(char_string2)-1),*) n ! last vertex

      nn = j
      nod(1) = m
      nod(2:j-1) = nod(2:j-1)
      nod(nn) = n

end subroutine


subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is a value between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is a value between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end

subroutine reverse_list(arrayIn,len)
!
! Given the list of few elements, and its size
! it reverses the order of elements and returns new 
! array.
!
  implicit none
  integer, intent(in) :: len
  integer, dimension(len), intent(inout) :: arrayIn

  integer,dimension(len) :: arrayOut

  integer :: i 

  do i=1,len
    arrayOut(len-i+1) = arrayIn(i)
  enddo

  arrayIn = arrayOut

end subroutine


subroutine find_opposite_face_or_pivot_edge(cellcon, faceNew, len, opposite, pivot)
!
! Oposite face to the current one in a hexahedral cell is the one
! which shares no vertex with a give face.
!
  implicit none
  integer :: i,j,l
  integer :: len
  integer, dimension(0:9) :: cellcon
  integer, dimension(len) :: faceNew
  logical :: opposite
  integer, dimension(2) :: pivot
  integer :: numOldVtx
  logical :: found_pivot

  numOldVtx = cellcon(0)

  ! Assume it is the opposite face and test against this assumption.
  opposite = .true.

  newface_loop: do i=1,len ! Take every vertex from the new face
    do j=1,numOldVtx ! and every vertex from the connectivity list

      if( faceNew(i) == cellcon(j) ) then
        ! We have found one sharing vertex - we call it first_pivot
        opposite = .false.
        pivot(1) = faceNew(i) 

        ! We should be able to find second_pivot vertex
        ! We take the next vertex in the new face list (index i+1) and check  
        ! if it is not shared with the old face. 
        ! Then we have found the pivot edge and of course second_pivot vertex.

        ! Assume we have found the new pivot at i+1 and test against this
        ! by checking that the next one is not shared node also:
        found_pivot = .true. 
        do l = 1,numOldVtx
          if( faceNew(i+1) == cellcon(l) ) found_pivot = .false.
        enddo

        if( found_pivot ) then

          if ( i+1 <= len ) then ! don't get out of index bounds for the faceNew array
            pivot(2) = faceNew(i+1)
          else
            pivot(2) = faceNew(1)
          endif
          exit newface_loop

        else

          cycle newface_loop

        endif

      endif

    enddo
  enddo newface_loop

end subroutine


subroutine cyclically_permute_nodes(cellcon, len, pivot)
implicit none

integer :: len
integer :: pivot
integer, dimension(len) :: cellcon
integer :: i
integer :: tmp
integer :: iter

! Case where there's nothing to cyclicaly_permute_nodes
! if( cellcon(1) == pivot ) return

iter = 0

do while( cellcon(1) /= pivot )

  iter = iter + 1

  if(iter >=len+1) then
    write(*,*) &
    "Something is wrong: Pivot not found in opposite face!"
    stop
  endif
  
  tmp = cellcon(1)

  do i=1,len-1
    cellcon(i) = cellcon(i+1)
  enddo

  cellcon(len) = tmp

enddo

return
end subroutine

integer function cell_type( nnodes )
implicit none
integer, intent(in) :: nnodes

if(nnodes == 8) then
  cell_type = 12 ! A hex in Paraview
elseif(nnodes == 6) then
  cell_type = 13 ! A prism in Paraview
elseif(nnodes == 5) then
  cell_type = 14 ! A pyramid in Paraview
elseif(nnodes == 4) then
  cell_type = 10 ! A tetra in Paraview
elseif(nnodes == 40) then
  cell_type = 9 ! A quad in Paraview
elseif(nnodes == 30) then
  cell_type = 5 ! A tetra in Paraview
else
  write(*,*) "Couldn't resolve cell type in cell_type function!"
  stop
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


end module utils