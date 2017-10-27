 program gmshToCappuccino
!
! Description:
!
! Converts Gmsh .msh files to one for use in Cappuccino solver.
!
! Author:
!   Nikola Mirkov
! 
! Date:
!   5. November 2015.
!

 use util
 use qsort_c_module
 implicit none

 integer, parameter :: dp = kind(1.0d0)
 integer, parameter :: largeInt = 999999999 ! safe, largest int32 (4-byte) integer is 2147483647
 integer, parameter :: nonoel = 8 ! no. of nodes in element-here Hex
 integer, parameter :: nonofa = 4 ! no, of nodes in element face-here Hex
 integer, parameter :: nofaelmax     = 6 ! no. of faces in element
 integer, parameter :: nofael_NTYPE4 = 4 ! no. of faces in element-here Tet
 integer, parameter :: nofael_NTYPE5 = 6 ! no. of faces in element-here Hex
 integer, parameter :: nofael_NTYPE6 = 5 ! no. of faces in element-here Prism
 integer, parameter :: nofael_NTYPE7 = 5 ! no. of faces in element-here Pyramid
 integer, parameter :: five = 5 ! Simply number five

 integer :: iarg
 integer :: iargc
 integer :: ios
 integer :: num_arg

 integer :: i,j,k,l
 integer :: iel, jel
 integer :: indx, indxEnd
 integer :: iface, jface
 integer :: nel    ! no. of elements in mesh
 integer :: ndim   ! dimension of the problem. 3 for 3D.
 integer :: nonome ! no. of nodes in mesh
 integer :: numhex ! No. of Hex elements 
 integer :: numpri ! No. of Prism elements
 integer :: numtet ! No. of Tet elements 
 integer :: numpyr ! No. of Pyr elements 
 integer :: nofalast ! No. of faces in the last element in the list
 integer :: njump
 integer :: nbc
 integer :: nInnerFaces
 integer :: nBndryFaces
 integer :: lenPN
 integer :: ivrtx

 character ( len = 255 ) prefix
 character ( len = 255 ) input_filename

 logical :: havePhysicalNames

! Gmsh related
 integer :: NUMNP, NELEM
 integer :: NP
 integer :: NE, NTYPE, NODE(8), NTAGS, TAGS(4)
 character(len=32) :: inLine
 character(len=10) :: numChar,numCharr
 real(dp), dimension(3) :: XCOO

 character(len=12), dimension(:), allocatable :: PhysicalName
 integer, dimension(:), allocatable :: PhysicalNameTag,PhysicalNameNum
 integer, dimension(:), allocatable :: cface_indx_start
 integer, dimension(:,:), allocatable :: faceVertices ! size[4,value_of(cface(6,nel))] or [4, 6*nel], a kolone su vrednosti za constituting_vertices eg.([ 1, 45, 23, 78]).
 integer, dimension(:,:), allocatable :: faceVerticesUnsorted
 integer, dimension(:,:), allocatable :: bFaceVertices
 integer, dimension(:,:), allocatable :: bFaceVerticesUnsorted
 integer, dimension(:), allocatable :: owner
 integer, dimension(:), allocatable :: neighbour



! 1) Intro
!+-----------------------------------------------------------------------------+
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  write ( *, '(a)' ) 'gmshToCappuccino'
  write ( *, '(a)' ) '  A preprocessor program for the Cappuccino code.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Reads mesh in Gmsh .msh format and prepares    '
  write ( *, '(a)' ) '  arrays encoding mesh sparsity pattern.         '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
!
!  Get the number of command line arguments.
!
  num_arg = iargc ( )

  if ( num_arg < 1 ) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'Enter the filename prefix:'
    read ( *, '(a)', iostat = ios ) prefix

    if ( ios /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Fatal error!'
      write ( *, '(a)' ) '  Unexpected read error!'
      stop
    end if

  else

    iarg = 1

    call getarg ( iarg, prefix )

  end if
!
!  Create the filenames.
!
  input_filename = trim ( prefix ) // '.msh'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Input file is "' // trim ( input_filename ) // '".'
  write ( *, '(a)' ) ' '
!+-----------------------------------------------------------------------------+



! > Read input from Gmsh mesh file
!+-----------------------------------------------------------------------------+
  open(unit=4,file=input_filename,status='old')
  rewind 4

  havePhysicalNames = .false.

! Skip header
! Aproach the line where $PhysicalNames are:
  do i=1,6
    read(4,'(a)') inLine
    if (inline(1:14).ne.'$PhysicalNames') then
      cycle
    else
      havePhysicalNames = .true.
      exit
    endif
  enddo

  if(.not.havePhysicalNames) then
    write(*,*) "File doesn't have $PhysicalNames defined. Quitting..."
    stop
  endif

! First thing we read is below line:
! k  - number of PhysicalNames
  read(4,*) k

  allocate( PhysicalName(k) )
  allocate( PhysicalNameNum(k) )
  allocate( PhysicalNameTag(k) )

  PhysicalNameNum = 0

  write( *, '(3x,a)' ) 'Physical Names:'

  lenPN = 0
  do i=1,k
    read(4,*) l,PhysicalNameTag(i),PhysicalName(i)
    write(*,'(3x,i2,1x,a)') PhysicalNameTag(i),PhysicalName(i)

    ! For boundary surfaces l=2, for inner domain l=3
    if (l.eq.2) lenPN = lenPN + 1
  enddo


  ! Skip $EndPhysicalNames and $Nodes lines
  do i=1,2
    read(4,*)
  enddo


! First thing we read is below line:
! NUMNP  - nuber of nodes in the mesh
  read(4,*) NUMNP

  write( *, '(a)' ) ' '
  write( *, '(3x,a)') 'NumPoints'
  write( *, '(1(1X,I9))') NUMNP
  write( *, '(a)' ) ' '


  nonome = NUMNP
  ndim = 3 !NDFCD
  if( ndim .ne. 3 ) then
      write(*,*) 'Fatal error: Only 3D meshes are accepted!'
      stop
  end if 


!open polyMesh format file: 'points'
  open(unit=7,file='points')
  rewind 7


!$Nodes
! Read nodal coordinates
  do i=1,nonome

    read(4,*) NP,(XCOO(k),k=1,ndim)

!...#Write to polyMesh format file: 'points'
    write(7,'(3E20.11)') (XCOO(k),k=1,ndim)

  enddo

!$EndNodes
! Skip rows '$EndNodes' and '$Elements'
  do i=1,2
    read(4,*)
  enddo

! NELEM  - nuber of cells in the mesh
  read(4,*) NELEM
  njump=0

!
! > OPEN polyMesh format file: 'cells',
!
  open(unit=12,file='cells')
  rewind 12
!
! > OPEN polyMesh format file: 'boundary',
!
  open(unit=8,file='boundary')
  rewind 8


! Read elements
  boundary_face_loop: do i=1,NELEM

!$Elements
! Format:
!$Elements
!number-of-elements
!elm-number elm-type number-of-tags < tag > … node-number-list
!…
!$EndElements
!
!Variable 	Description
!NE 	Global element number (not required to be sequential or continuous)
!NTYPE 	Element geometry type:
!	1 = 2-node line.
!	2 = 3-node triangle.
!	3 = 4-node quadrangle.
!	4 = 4-node tetrahedron.
!	5 = 8-node hexahedron.
!	6 = 6-node prism.
!	7 = 5-node pyramid.
!       ... (Whole list and descritpion is in README section: gmsh-msh-format.)
!NTAGS 	Number of tags for the element
!TAGS   List of tags for the element
!NODE 	List of nodes that define the element

  read(4,*) NE, NTYPE, NTAGS, (TAGS(k), k=1,NTAGS), (NODE(k), k=1,noel(NTYPE))


! > Maybe we have face describing boundary conditions...
  if (NTYPE.eq.2 .or. NTYPE.eq.3) then ! Triangular or Quad boundary face:
    njump = njump+1
    write(8,*) TAGS(2),NTYPE,(NODE(k), k=1,noel(NTYPE))
    call addTaggedBFaceToTotalNumber( TAGS(1), PhysicalNameTag, PhysicalNameNum, lenPN )
  elseif(NTYPE.eq.1) then
    print*, 'We have hit an edge in the elements list, we need faces or cells!'
    stop
  else
    if (i.eq.1) njump = 0 
    exit boundary_face_loop
  endif

end do boundary_face_loop

 ! Return one line back
 backspace(4)

  nel = NELEM-njump


! > We can read cell now

  write( *, '(3x,a)') 'NumCells'
  write( *, '(1(1X,I9))') nel
  write( *, '(a)' ) ' '



  allocate ( faceVertices(nonofa+1,nofaelmax*nel) )
  allocate ( cface_indx_start(nel+1) )

! Initialize
  faceVertices(:,:) = 0


  numhex = 0
  numpri = 0
  numtet = 0
  numpyr = 0

  cface_indx_start(1) = 1

! Read elements
  element_loop: do i=1,nel

!$Elements
! Format:
!$Elements
!number-of-elements
!elm-number elm-type number-of-tags < tag > … node-number-list
!…
!$EndElements
!
!Variable 	Description
!NE 	Global element number (not required to be sequential or continuous)
!NTYPE 	Element geometry type:
!	1 = 2-node line.
!	2 = 3-node triangle.
!	3 = 4-node quadrangle.
!	4 = 4-node tetrahedron.
!	5 = 8-node hexahedron.
!	6 = 6-node prism.
!	7 = 5-node pyramid.
!       ... (Whole list and descritpion is in README section: gmsh-mah-format.)
!NTAGS 	Number of tags for the element
!TAGS   List of tags for the element
!NODE 	List of nodes that define the element

  read(4,*) NE, NTYPE, NTAGS, (TAGS(k), k=1,NTAGS), (NODE(k), k=1,noel(NTYPE))

  cface_indx_start(i+1) = cface_indx_start(i) + nofael(NTYPE)

!
! > Write into 'cells' polyMesh file
!
  write(12,'(I1,1X,4I8:(4I8:))') ntypeCappuccino(NTYPE),(NODE(k), k=1,noel(NTYPE))


 if (NTYPE.eq.4 .or. NTYPE.eq.11 .or. NTYPE.eq.29 ) then
! ELEMENT is 4-node TETRAHEDRON
! Tetrahedron:                          Tetrahedron10:
! 
!                    v
!                  .
!                ,/
!               /
!            2                                     2                              
!          ,/|`\                                 ,/|`\                          
!        ,/  |  `\                             ,/  |  `\       
!      ,/    '.   `\                         ,6    '.   `5     
!    ,/       |     `\                     ,/       8     `\   
!  ,/         |       `\                 ,/         |       `\ 
! 0-----------'.--------1 --> u         0--------4--'.--------1
!  `\.         |      ,/                 `\.         |      ,/ 
!     `\.      |    ,/                      `\.      |    ,9   
!        `\.   '. ,/                           `7.   '. ,/     
!           `\. |/                                `\. |/       
!              `3                                    `3        
!                 `\.
!                    ` w
 

   ! How many faces before
   jface  = numtet*nofael_NTYPE4 + numhex*nofael_NTYPE5 + numpri*nofael_NTYPE6 + numpyr*nofael_NTYPE7 

   do iface = 1,nofael_NTYPE4
     faceVertices(:, jface+iface ) = get_face_vertices_NTYPE4(iface,NODE,noel(NTYPE))
     faceVertices(5,jface+iface) = i ! = element number
   enddo

   numtet = numtet + 1

 elseif (NTYPE.eq.5 .or. NTYPE.eq.17 .or. NTYPE.eq.12) then
! ELEMENT is 8-node HEX
! Hexahedron:             Hexahedron20:          Hexahedron27:
!
!        v
! 3----------2            3----13----2           3----13----2     
! |\     ^   |\           |\         |\          |\         |\    
! | \    |   | \          | 15       | 14        |15    24  | 14  
! |  \   |   |  \         9  \       11 \        9  \ 20    11 \  
! |   7------+---6        |   7----19+---6       |   7----19+---6 
! |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23| 
! 0---+---\--1   |        0---+-8----1   |       0---+-8----1   | 
!  \  |    \  \  |         \  17      \  18       \ 17    25 \  18
!   \ |     \  \ |         10 |        12|        10 |  21    12| 
!    \|      w  \|           \|         \|          \|         \| 
!     4----------5            4----16----5           4----16----5 
!


   ! How many faces before
   jface  = numtet*nofael_NTYPE4 + numhex*nofael_NTYPE5 + numpri*nofael_NTYPE6 + numpyr*nofael_NTYPE7

   do iface = 1,nofael_NTYPE5
     faceVertices(:, jface+iface ) = get_face_vertices_NTYPE5(iface,NODE,noel(NTYPE))
     faceVertices(5,jface+iface) = i ! = element number
   enddo

   numhex = numhex + 1

 elseif (NTYPE.eq.6 .or. NTYPE.eq.18 .or. NTYPE.eq.13) then
! ELEMENT is 6-node PRISM
!  Prism:                      Prism15:               Prism18:
!
!            w
!            ^
!            |
!            3                       3                      3        
!          ,/|`\                   ,/|`\                  ,/|`\      
!        ,/  |  `\               12  |  13              12  |  13    
!      ,/    |    `\           ,/    |    `\          ,/    |    `\  
!     4------+------5         4------14-----5        4------14-----5 
!     |      |      |         |      8      |        |      8      | 
!     |    ,/|`\    |         |      |      |        |    ,/|`\    | 
!     |  ,/  |  `\  |         |      |      |        |  15  |  16  | 
!     |,/    |    `\|         |      |      |        |,/    |    `\| 
!    ,|      |      `\        10     |      11       10-----17-----11
!  ,/ |      0      | `\      |      0      |        |      0      | 
! u   |    ,/ `\    |    v    |    ,/ `\    |        |    ,/ `\    | 
!     |  ,/     `\  |         |  ,6     `7  |        |  ,6     `7  | 
!     |,/         `\|         |,/         `\|        |,/         `\| 
!     1-------------2         1------9------2        1------9------2 
!

   ! How many faces before
   jface  = numtet*nofael_NTYPE4 + numhex*nofael_NTYPE5 + numpri*nofael_NTYPE6 + numpyr*nofael_NTYPE7

   do iface = 1,nofael_NTYPE6
     faceVertices(:, jface+iface ) = get_face_vertices_NTYPE6(iface,NODE,noel(NTYPE))
     faceVertices(5,jface+iface) = i ! = element number
   enddo

   numpri = numpri + 1


 elseif (NTYPE.eq.7 .or. NTYPE.eq.19 .or. NTYPE.eq.14) then
! ELEMENT is 5-node PYRAMID
!  Pyramid:                     Pyramid13:                   Pyramid14:
! 
!                4                            4                            4
!              ,/|\                         ,/|\                         ,/|\
!            ,/ .'|\                      ,/ .'|\                      ,/ .'|\
!          ,/   | | \                   ,/   | | \                   ,/   | | \
!        ,/    .' | `.                ,/    .' | `.                ,/    .' | `.
!      ,/      |  '.  \             ,7      |  12  \             ,7      |  12  \
!    ,/       .' w |   \          ,/       .'   |   \          ,/       .'   |   \
!  ,/         |  ^ |    \       ,/         9    |    11      ,/         9    |    11
! 0----------.'--|-3    `.     0--------6-.'----3    `.     0--------6-.'----3    `.
!  `\        |   |  `\    \      `\        |      `\    \     `\        |      `\    \
!    `\     .'   +----`\ - \ -> v  `5     .'        10   \      `5     .' 13     10   \
!      `\   |    `\     `\  \        `\   |           `\  \       `\   |           `\  \ 
!        `\.'      `\     `\`          `\.'             `\`         `\.'             `\` 
!           1----------------2            1--------8-------2           1--------8-------2
!                     `\
!                        u
!

   ! How many faces before
   jface  = numtet*nofael_NTYPE4 + numhex*nofael_NTYPE5 + numpri*nofael_NTYPE6 + numpyr*nofael_NTYPE7

   do iface = 1,nofael_NTYPE7
     faceVertices(:, jface+iface ) = get_face_vertices_NTYPE7(iface,NODE,noel(NTYPE))
     faceVertices(5,jface+iface) = i ! = element number
   enddo

   numpyr = numpyr + 1

 elseif (NTYPE.eq.2) then
    cycle ! We have hit the triangular boundary face.
 elseif (NTYPE.eq.3) then
     cycle ! We have hit the quad boundary face.
 else
      write(*,*) 'Fatal error: Non-existing cell type!'
      stop
 endif

 if (i.eq.nel) then
   nofalast = nofael(NTYPE)
 endif

end do element_loop

!
! > CLOSE polyMesh format file: 'cells',
!
  close(12)



!
! Report
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Gmsh file:'
  write ( *, '(a,I7)' ) '  Total no. of HEX cells: ', numhex
  write ( *, '(a,I7)' ) '  Total no. of TET cells: ', numtet
  write ( *, '(a,I7)' ) '  Total no. of PYR cells: ', numpyr
  write ( *, '(a,I7)' ) '  Total no. of PRI cells: ', numpri
  write ( *, '(a)' )    ' +--------------------------------='
  write ( *, '(a,I7)' ) '  Total no. of cells: ', numhex+numtet+numpyr+numpri
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Normal end of reading .msh file.'
  write ( *, '(a)' ) ' '


!+-----------------------------------------------------------------------------+
! > END: Read input from Gmsh mesh file.


  allocate ( faceVerticesUnsorted(nonofa,nofaelmax*nel) )
  faceVerticesUnsorted = faceVertices(1:nonofa, :)

  write ( *, '(a)' ) 'Sort faceVertices: Begin'
    do i = 1,nel*nofaelmax
      call sortIntArray(faceVertices(1:nonofa,i),nonofa)
    enddo
  write ( *, '(a)' ) 'Sort faceVertices: End'
  write ( *, '(a)' ) ' '

  write ( *, '(a)' ) 'First 20 lines of faceVertices array: '
  write ( *, '(a)' ) ' '
    do i = 1,20
      write ( *, * ) faceVertices(:,i)
    enddo
  write ( *, '(a)' ) ' '

!
! Boundary conditions: Continue
!

!+-----------------------------------------------------------------------------+

!
! > Rewind polyMesh format file: 'boundary',
!
  rewind 8

 ! Odavde ponovo procitaj 'boundary' nadji koji celiji pripada i lepo zapisi u isti fajl 
 nbc = njump

 allocate ( bFaceVertices(nonofa,nbc) )
 allocate ( bFaceVerticesUnsorted(nonofa,nbc) )

 bFaceVertices(:,:) = 0
 bFaceVerticesUnsorted(:,:) = 0

! Read boundary 
  do i=1,nBC
    read(8,*) TAGS(2),NTYPE,(bFaceVertices(k,i), k=1,noel(NTYPE))
    ! Back-up
    bFaceVerticesUnsorted(:,i) = bFaceVertices(:,i)
    call sortIntArray(bFaceVertices(:,i),nonofa)
  end do 

!
! > Rewind polyMesh format file: 'boundary',
!
  rewind 8
  
  write ( *, '(a)' ) ' Find boundary cells:'
!1=============================================================================1!
  nBndryFaces = 0

  bc_loop: &
  do jface = 1,nbc

  ! Progress status:
  ! write(*,'(11a,f5.1,a)',advance='no') char(8),char(8),char(8),char(8),char(8),char(8),char(8),&
  ! char(8),char(8),char(8),char(8),float(jface)/float(nbc)*100,'% done'
 
  do iel = 1,nel

    indx = cface_indx_start(iel)
    indxEnd = cface_indx_start(iel+1)-1


    do iface = indx, indxEnd !----------------------------------------------->

    ivrtx = 1
    if ( faceVertices(1,iface) == 0 ) ivrtx = 2

      if ( faceVertices(ivrtx,iface) .eq. bFaceVertices(ivrtx,jface) ) then 
      ! znaci dele zajednicki vertex - interesantno idemo dalje...

        if ( faceVertices(ivrtx+1,iface) .eq. bFaceVertices(ivrtx+1,jface) ) then 
        ! znaci dele zajednicku stranicu, jos bolje! Idemo dalje...

          if ( faceVertices(ivrtx+2,iface) .eq. bFaceVertices(ivrtx+2,jface) ) then 
          ! To je to, imaju tri zajednica vertexa znaci da su susedi zapisujemo to:

            ivrtx = 1
            if ( bFaceVertices(4,jface) == 0 ) ivrtx = 2
            write(8,*) faceVertices(5,iface), (five-ivrtx), bFaceVerticesUnsorted(1:nonofa-ivrtx+1,jface)

            nBndryFaces = nBndryFaces + 1

           cycle bc_loop

          endif
        endif
      endif

    enddo  !<--------------------------------------------------------------------

  enddo 

  enddo bc_loop
!1=============================================================================1!
  write ( *, '(a)' ) ' Done.'
  write ( *, '(a)' ) ' '

! We don't need this anymore...
  deallocate ( bFaceVerticesUnsorted )
  deallocate ( bFaceVertices )


!+-----------------------------------------------------------------------------+

!
! > OPEN polyMesh format file: 'faces', 'owner', 'neighbour'.
!
  open(unit=9,file='faces')
  open(unit=10,file='owner')
  open(unit=11,file='neighbour')


  rewind 9
  rewind 10
  rewind 11

! > QUICKSORT over vertices:

  write ( *, '(a)' ) ' Quicksort faceVertices: Begin'

  ! Do quicksort on a first vertex  - should speed things up considerably
    call QsortC(faceVertices(2,:),faceVertices(3,:),faceVertices(4,:),faceVertices(1,:),faceVertices(5,:), &
                faceVerticesUnsorted(1,:),faceVerticesUnsorted(2,:),faceVerticesUnsorted(3,:),faceVerticesUnsorted(4,:))

   i=1
   k=i+1
   first_index_loop:do iel=1,nofaelmax*nel
   if(faceVertices(2,i)==faceVertices(2,k)) then
    k=k+1
   else
    ! Ready to Go:
    call QsortC(faceVertices(3,i:k-1),faceVertices(4,i:k-1),faceVertices(1,i:k-1),faceVertices(2,i:k-1),faceVertices(5,i:k-1), &
      faceVerticesUnsorted(1,i:k-1),faceVerticesUnsorted(2,i:k-1),faceVerticesUnsorted(3,i:k-1),faceVerticesUnsorted(4,i:k-1))
    i=k
    k=i+1
    cycle first_index_loop
   endif
   enddo first_index_loop

   i=1
   k=i+1
   second_index_loop:do iel=1,nofaelmax*nel
   if(faceVertices(3,i)==faceVertices(3,k)) then
    k=k+1
   else
    ! Ready to Go:
    call QsortC(faceVertices(4,i:k-1),faceVertices(1,i:k-1),faceVertices(2,i:k-1),faceVertices(3,i:k-1),faceVertices(5,i:k-1), &
      faceVerticesUnsorted(1,i:k-1),faceVerticesUnsorted(2,i:k-1),faceVerticesUnsorted(3,i:k-1),faceVerticesUnsorted(4,i:k-1))
    i=k
    k=i+1
    cycle second_index_loop
   endif
   enddo second_index_loop

! > END: QUICKSORT over vertices.

  write ( *, '(a)' ) ' Quicksort faceVertices: End'
!
! > Main loops for defining sparsity pattern: iel_loop, jel_loop
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Find cells with common inner face:'
  

!2=============================================================================2!

  iel = 1 
  nInnerFaces = 0

  indxEnd = size(faceVertices,2)-1
  iface_loop: do iface=1,indxEnd

    ! Progres status:
    ! write(*,'(11a,f5.1,a)',advance='no') char(8),char(8),char(8),char(8),char(8),char(8),char(8),&
    ! char(8),char(8),char(8),char(8),float(iface)/float(indxEnd)*100,'% done'

    if ( faceVertices(2,iface) == 0 ) cycle iface_loop

    j = 2
    if (iface .eq. indxEnd) j=1
    do jface = iface+1, iface+j

      ivrtx = 1
      if ( faceVertices(1,iface) == 0 ) ivrtx = 2 ! tj. kad je face trougao.

      if ( faceVertices(2,iface) .eq. faceVertices(2,jface) ) then 
      ! znaci dele zajednicki vertex - interesantno idemo dalje...

        if ( faceVertices(3,iface) .eq. faceVertices(3,jface) ) then
        ! znaci dele zajednicku ivicu, jos bolje! Idemo dalje...

          if ( faceVertices(4,iface) .eq. faceVertices(4,jface) ) then
          ! E to je to, imaju tri zajednica vertexa znaci da su susedi zapisujemo to

            iel = faceVertices(5,iface)
            jel = faceVertices(5,jface)

              ! write face into the face file:
              write(9,'(4i8)') faceVerticesUnsorted(1:nonofa-ivrtx+1,iface)

              if (iel < jel) then
                ! write owner:
                write(10,'(i8)') iel
                ! write neighbour
                write(11,'(i8)') jel 
              else
                ! write owner:
                write(10,'(i8)') jel
                ! write neighbour
                write(11,'(i8)') iel 
              endif

              nInnerFaces = nInnerFaces + 1
        
              cycle iface_loop

          endif
        endif
      endif
    enddo
  enddo iface_loop

 ! We don't need this anymore...
 deallocate ( faceVertices )
 deallocate ( faceVerticesUnsorted )
 deallocate ( cface_indx_start )

!2=============================================================================2!
 
  write ( *, '(a)' ) ' Done.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) ' Number of Inner faces = ', nInnerFaces

!
!  > Rewind polyMesh format file: 'faces', 'owner', 'neighbour'.
!
  rewind 9
  rewind 10
  rewind 11

!
! > Allocate arrays for owner, neighbour ad face vertices
!
  allocate ( owner( nInnerFaces) )
  allocate ( neighbour( nInnerFaces) )
  allocate ( faceVertices( nonofa, nInnerFaces) )


  do iface = 1,nInnerFaces
    read(9,'(4i8)') faceVertices(:,iface)
    read(10,'(i8)') owner(iface) 
    read(11,'(i8)') neighbour(iface)
  end do

!
! > QUICKSORT the arrays
!

! Do quicksort on owner array
  call QsortC2( owner(:),neighbour(:), &
                faceVertices(1,:),faceVertices(2,:),faceVertices(3,:),faceVertices(4,:) )

 i=1
 k=i+1
 index_loop:do iel=1,nInnerFaces
 if(owner(i)==owner(k)) then
  k=k+1
 else
  ! Ready to Go:
    call QsortC2( neighbour(i:k-1),owner(i:k-1), &
                  faceVertices(1,i:k-1),faceVertices(2,i:k-1),faceVertices(3,i:k-1),faceVertices(4,i:k-1) ) 
  i=k
  k=i+1
  cycle index_loop
 endif
 enddo index_loop

!
! > Now write polyMesh files again with sorted arrays
!
  rewind 9
  rewind 10
  rewind 11

  ! Write header information
  call i4_to_s_left ( nonome, numChar ) 
  write(10,'(a,a)') trim(numChar),' points;'

  call i4_to_s_left ( nel, numChar ) 
  write(10,'(a,a)') trim(numChar),' cells;'

  call i4_to_s_left ( (nInnerFaces + nBndryFaces), numChar ) 
  write(10,'(a,a)') trim(numChar),' faces;'

  call i4_to_s_left ( nInnerFaces, numChar ) 
  write(10,'(a,a)') trim(numChar),' internal faces.'
  
  do iface = 1,nInnerFaces

    ivrtx = 1 ! For a quad face ...
    if ( faceVertices(4,iface) == 0 ) ivrtx = 2 !..if it is actually a tri face.

    call i4_to_s_left ( (five-ivrtx), numChar ) 
    write(9,'(a,4i8)') trim(numChar), faceVertices(1:nonofa-ivrtx+1,iface)

    call i4_to_s_left ( owner(iface), numChar ) 
    write(10,'(a)') numChar

    call i4_to_s_left ( neighbour(iface), numChar ) 
    write(11,'(a)') numChar

  end do



! !
! ! > Quicksort the boundary cells in ascending order
! !

!   rewind 8

!   allocate ( bFaceVertices( nonofa, nBndryFaces) )

!   do i=1,nBndryFaces
!     read(8,*) owner(i), ivrtx, bFaceVertices(:,i)
!   enddo

!   rewind 8

! !
! ! > Quicksort the boundary cells in ascending order
! !
!   write ( *, '(a)' ) ' '
!   write ( *, '(a)' ) ' Quicksort bfaceVertices: Begin'
!   call QsortC3( owner(1:nBndryFaces), &
!                 bFaceVertices(1,:),bFaceVertices(2,:),bFaceVertices(3,:),bFaceVertices(4,:) )
!   write ( *, '(a)' ) ' Quicksort bfaceVertices: End'
!   write ( *, '(a)' ) ' '

! !
! ! > Write to 'boundary' file again
! !
!   do i=1,nBndryFaces
!     ivrtx = 1
!     if ( bFaceVertices(4,i) == 0 ) ivrtx = 2
!     write(8,*) owner(i), (five-ivrtx), bFaceVertices(1:nonofa-ivrtx+1,i)
!   enddo

!   deallocate ( bFaceVertices )

  deallocate ( owner )
  deallocate ( neighbour )
  deallocate ( faceVertices )


!
! > Write boundary faces to 'owner' and 'faces' files
!
  rewind 8
  do i=1,nBndryFaces
    read(8,*) iel, k, node(1:k)
    call i4_to_s_left ( iel, numChar ) 
    write(10,'(a)') numChar
    call i4_to_s_left ( k, numChar ) 
    write(9,'(a,4i8)') trim(numChar), node(1:k)
  enddo

!
! > Write boundary conditions based on Physical names
!
  rewind 8

  write(8,'(a)') '#bctype nFaces startFace'

  l = nInnerFaces 
  do i=1,lenPN
    call i4_to_s_left ( l, numChar )
    call i4_to_s_left ( PhysicalNameNum(i), numCharr )
    write(8,'(a,1x,a,1x,a)') trim( PhysicalName(i) ), trim( numCharr ), numChar
    l = l + PhysicalNameNum(i)
  enddo

!
!  > CLOSE polyMesh format file: 'points','boundary' faces', 'owner', 'neighbour'.
!
  close (7)
  close (8)
  close (9)
  close (10)
  close (11)


!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Preprocessor:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

 end program gmshToCappuccino