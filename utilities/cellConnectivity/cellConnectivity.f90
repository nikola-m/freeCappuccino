!*******************************************************************************
program cellConnectivity
!
!  Description:
!    Part of the writer for VTK unstructured .vtu format.
!
!  Discussion:
!    Main difficulty is in creating cell connectivity data from face
!    based data structure of polyMesh format.
!    This utility should be used if one is provided with mesh in OpenFOAM 
!    polyMesh format. In other cases, such as when mesh data is given in 
!    Gmsh .msh format, we already have cell connectivity and preparing data 
!    for postprocessing is much easier.
!
!  Date:
!    15/01/2018
!
!  Author:
!    Nikola Mirkov largeddysimulation@gmail.com or nmirkov@vin.bg.ac.rs
!
!*******************************************************************************
!
  use utils

  implicit none


  ! Logical defining are we reading native Cappuccino file format or OpenFOAM polymesh format.
  logical :: native_mesh_files

  ! Mesh file units
  integer :: faces_file, owner_file, neighbour_file, cells_file

  character(len=80) :: line_string
  character(len=15) :: char_string
  character(len=1) :: ch

  integer, parameter :: nomax = 8 ! Change this if necessary, for size.
  integer :: iface
  integer :: i,j,k,l
  integer :: len,len2
  integer :: inp,inn,in
  integer :: nnodes
  integer :: ctype
  integer :: numCells, numFaces, numInnerFaces, numBoundaryFaces, numBFace
  integer :: node(nomax) 
  logical :: opposite
  integer :: pivot(2),pivt
  integer, dimension(:,:), allocatable :: connectivity

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
  k=0
  l=0
  do i=1,17
    if (i==13) then
      read(owner_file,*) char_string,line_string
      do j=10,len_trim(line_string)
        if (line_string(j:j+5)=='nCells') then
          k=j+7
        endif
        if (line_string(j:j+5)=='nFaces') then
          l=j-2
        endif
      end do
      read(line_string(k:l),*) numCells
    else
      read(owner_file,*) ch
    endif
  end do
  read(owner_file,*) numFaces
  read(owner_file,*) ch ! reads "("


  ! write(*,*) "numCells: ", numCells
  ! write(*,*) "numFaces: ", numFaces

  ! 'neighbour' file
  do i=1,17
    read(neighbour_file,*) ch
  end do
  read(neighbour_file,*) numInnerFaces
  read(neighbour_file,*) ch ! reads "("


  numBoundaryFaces = numFaces - numInnerFaces

  ! 'faces' file
  do i=1,18
    read(faces_file,*) ch
  end do

  ! write(*,*) "numInnerFaces: ", numInnerFaces

  ! Allocate array with cell connectivity information
  allocate ( connectivity(0:9, 0:numCells+numBoundaryFaces-1) )

  ! Initialize connectivity
  connectivity  = -1

  ! Test important routine:
  ! node(1:4) = ([1,2,3,4])
  ! call cyclically_permute_nodes( node(1:4), 4, 4 )
  ! write(*,*) node(1:4)
  ! stop

  ! Loop over faces_file
  face_loop: do iface=1,numFaces

    read(owner_file,*) inp

    if (iface <= numInnerFaces ) then
      read(neighbour_file,*) inn
    endif

    node = 0

    ! Read line in 'faces' file
    if (native_mesh_files) then
      read( faces_file, * ) nnodes, (node(k), k=1,nnodes)
    else ! OpenFOAM polyMesh
      call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)
    endif


    ! Write boundary faces at the end of connectivity array
    if (iface > numInnerFaces ) then

      numBFace = numCells + (iface - numInnerFaces) - 1
       
      ! nnodes*10 - That way I will get 30 for tri and 40 for quad. 
      ! Otherwise I would have collision of qaud and tet bacuse they have the same number of nodes.
      connectivity(0, numBFace) = cell_type( 10*nnodes ) 

      connectivity(1:nnodes, numBFace) = node(1:nnodes)

    endif


    sides_loop: do l=1,2

      ! Repeat the same proces for both owner 'inp'  and neighbour 'inn', only that we have
      ! neighbours for iface =< numInnerFaces.
      if ( l == 1 ) then
        in = inp
        ! Write owner cell nodes but reverse the order of nodes in face. 
        call reverse_list(node(1:nnodes), nnodes) 

      elseif( iface <= numInnerFaces ) then
        in = inn

      else
        cycle face_loop

      endif

      if ( connectivity(0,in) == -1 ) then

       ! Nothing is yet written for this cell

        connectivity(0,in) = nnodes
   
        ! Write owner cell nodes but reverse the order of nodes in face.
        ! if ( l == 1 ) call reverse_list(node(1:nnodes), nnodes)


        connectivity(1:nnodes,in) = node(1:nnodes)


      elseif ( connectivity(9,in) == 999 .and. connectivity(9,inn) == 999) then
 
        cycle face_loop

      elseif ( connectivity(9,in) == 999 .and. connectivity(9,inn) /= 999) then
 
        cycle sides_loop

      else

        ! We have base face in place so we should try to find the opposite face suspecting that we are dealing with
        ! hex or prism cell, or, if not that, to find pivot edge.

        call find_opposite_face_or_pivot_edge(connectivity(:,in), node(1:nnodes), nnodes, opposite, pivot)

        ! if(.not.opposite) call find_pivot_edge(connectivity(:,in), node(1:nnodes), nnodes, pivot)

        if(opposite) then 

          ! 1)
          ! We have found opposite face to the base face, write the face in the connectivity list.
          ! But before check do we already have the pivot in place.

          if( connectivity(9,in) == 100 ) then

          ! 1.1) We have pivot in place, grab it
            pivt = connectivity( len+1 ,in)

            ! First cell bug thing
            if ( in == 0 ) call reverse_list(node(1:nnodes), nnodes)

            ! Use second pivot 'pivt' to cyc
            call cyclically_permute_nodes( node(1:nnodes), nnodes, pivt )

            len = connectivity(0,in)
            connectivity(len+1:len+nnodes,in) = node(1:nnodes)

            ! Write nuber of nodes in cell at place with index 0, which will be used to define cell type, 12 for HEX, 10 for TET, etc..
            connectivity( 0,in ) = len+nnodes

            ! Seal the cell - we did it 
            connectivity( 9,in ) = 999
          

          else

            ! 1.2) We don't have pivot edge on place
            ! We still have to cyclicaly permute nodes later...

            ! First cell bug thing
            if ( in == 0 ) call reverse_list(node(1:nnodes), nnodes)
            
            len = connectivity( 0,in )

            connectivity( len+1:len+nnodes, in ) = node( 1:nnodes )

          endif

        elseif ( .not.opposite .and. connectivity(9,in) == 100 ) then

          ! 2)
          ! We've got a face that is not the opposite face but we already have second pivot on place
          ! and we don't want to waste our computational resources so we just cycle to next face in line 
          if (l == 1) then
            cycle sides_loop
          else
            cycle face_loop
          endif

        else

          ! 3)
          ! We haven't got the opposite face, but, we're going to set pivot edge. We are going to leave the sign 
          ! that we've done that so we don't have to repeat it in the future, by setting:
          connectivity(9,in) = 100

          ! Use the first pivot to cyclically permute nodes in the base face
          len = connectivity(0,in)
          call cyclically_permute_nodes(connectivity(1:len, in), len, pivot(1))


          ! Use second pivot to permute the rest if we allready have 
          ! opposite face on its place,
          ! or,
          ! just write the second pivot after the base face nodes in the connectivity list.
   
          ! Test - if we have allready found the opposite face, whether we are dealing with the hex or prism
          ! (the two cases where we need the oposite face), there will be something other than -1 in connectivity
          ! list index 6:
          if (connectivity(6,in) /= -1) then 

            ! 3.1)
            ! we have two faces, base and opposite, in place
      
            ! Length of the opposite face, initialize:
            len2 = 0
            ! Calculate incrementaly:
            do i=1,(8-len) ! cell nodes are only in places 1:8, 0 and 9 are reserved for other things, see above.
              if ( connectivity( len+i, in ) /= -1) then
                len2 = len2 + 1
              else
                exit
              endif
            enddo

            call cyclically_permute_nodes( connectivity( len+1:len+len2, in ), len2, pivot(2) )

            ! Write nuber of nodes in cell at place with index 0, which will be used to define cell type, 12 for HEX, 10 for TET, etc..
            connectivity( 0,in ) = len+len2

            ! Seal the cell - we did it 
            connectivity( 9,in ) = 999
 
          else 

            ! 3.2)
            ! we don't have the opposite face yet, just write down, second pivot vertex

            connectivity( len+1, in ) = pivot(2) 
   
          endif

        endif ! if opposite conditional

      endif ! conditional to: write base face, write oposite face, or get pivot edge and either cyclically permute base and 
            ! opposite face, if both in place,or just pemute base face and write down second pivot for future reference. 

    enddo sides_loop

  enddo face_loop


  do i=1,numCells

    nnodes = connectivity(0,i-1)

    write(cells_file,'(i2,1x,8i8)') cell_type( nnodes ), (connectivity( k, i-1 ), k=1,nnodes)

  enddo

  do i=numCells, numCells+numBoundaryFaces-1

    ctype = connectivity(0,i)

    write(cells_file,'(i2,1x,8i8)') ctype, (connectivity( k, i ), k=1,noel(ctype))

  enddo  

 close(neighbour_file)
 close(owner_file)
 close(faces_file)
 close(cells_file)

end program