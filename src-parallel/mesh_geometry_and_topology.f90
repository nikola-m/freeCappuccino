!
! Module for geometry definition on unstructured meshes.
!
module geometry

use types
use parameters, only: myid
use utils, only: get_unit, file_row_count, r8vec_print_some, i4vec_print2, i4_to_s_left

implicit none

! General mesh data
integer :: numNodes                           ! no. of nodes in mesh
integer :: numCells                           ! no. of cells in the mesh
integer :: numPCells                          ! no. of cells in the mesh + buffer cells
integer :: numFaces                           ! no. of INNER+BOUNDARY faces in the mesh
integer :: numInnerFaces                      ! no. of INNER cells faces in the mesh
integer :: numBoundaryFaces                   ! self explanatory
integer :: numTotal                           ! number of volume field values + number of boundary field values numCells+numBoundaryFaces
integer :: nnz                                ! no. of nonzeros in sparse system matrix

integer :: ninl                               ! No. of inlet boundary faces
integer :: nout                               ! No. of outlet boundary faces
integer :: nsym                               ! No. of symmetry boundary faces
integer :: nwal                               ! No. of wall boundary faces
integer :: npru                               ! No. of pressure outlet boundary faces
integer :: noc                                ! No. of O-C- domain cut boundary faces
integer :: ncyc                               ! No. of cyclic boundary faces
integer :: npro                               ! No. of processor boundary faces
integer :: nwali,nwala,nwalf                  ! No. of isothermal,adiabatic,heatflux boundary faces (for temperature eq.)
integer :: iOCBoundaries                      ! No. of O-C- boundaries, they come in pair
integer :: iCycBoundaries                     ! No. of cyclic boundaries, they come in pair

! In variable arrays, field variables defined at cell centres are written in positions from 1 to numCells, after that we write
! variable values for boundary faces. We need to know where are values defined at boundary faces in these arrays.
! That is the purpose of variables defined below.
! We have for example that inlet values for u-velocity component can be found at u( iInletStart+1 : iInletStart+ninl ), 
! Temperature at wall for aexample can be found at positions t( iWallStart+1 : iWallStart+nwal ), and so on...

integer :: iInletStart                        
integer :: iOutletStart                       
integer :: iSymmetryStart                     
integer :: iWallStart                        
integer :: iPressOutletStart                  
integer :: iOCStart                           
integer :: iCycStart
integer :: iProcStart                          

! Where the faces pertinet to certain boundary types are located:

integer :: iInletFacesStart                   
integer :: iOutletFacesStart                  
integer :: iSymmetryFacesStart                
integer :: iWallFacesStart                    
integer :: iPressOutletFacesStart             
integer :: iOCFacesStart                      
integer :: iCycFacesStart    
integer :: iProcFacesStart                 

integer, parameter :: nomax = 24              ! Max no. of nodes in face - determines size of some arrays, just change this if necessary.
real(dp), parameter :: tiny = 1e-30

! Logical defining are we reading native Cappuccino file format or OpenFOAM polymesh format.
logical :: native_mesh_files

! Mesh file units
integer :: points_file, faces_file, owner_file, neighbour_file, boundary_file, process_file

! Global number of nodes, cells, faces, inner faces, and boundary faces when summed from all domains
integer :: gloNodes,gloCells, gloFaces, gloIFaces, gloBFaces


! Mesh geometry

! Geometry parameters defined for mesh nodes [1:numNodes]
real(dp), dimension(:), allocatable :: x ,y ,z  ! Coordinates of mesh nodes 

! Geometry parameters defined cellwise [1:numCells]:
real(dp), dimension(:), allocatable :: xc,yc,zc     ! Coordinates of cell centers
real(dp), dimension(:), allocatable :: vol          ! Cell volume
real(dp), dimension(:), allocatable :: wallDistance ! Distance to the nearest wall - needed in some turb. models

! Geometry parameters defined for all (inner+boundary) cell faces [1:numFaces]
real(dp), dimension(:), allocatable :: arx, ary, arz   ! Cell face area x-, y- and z-component 
real(dp), dimension(:), allocatable :: xf, yf, zf      ! Coordinates of cell face center

! Geometry parameters defined for all inner cell faces [1:numInnerFaces]
real(dp), dimension(:), allocatable :: xpp, ypp, zpp   ! Coordinates of auxilliary points - owner cell 
real(dp), dimension(:), allocatable :: xnp, ynp, znp   ! Coordinates of auxilliary points - neighbour cell 
real(dp), dimension(:), allocatable :: facint          ! Interpolation factor 
!real(dp), dimension(:), allocatable :: dpn            ! Distance between neigbor cell centers [1:numInnerFaces]

! Geometry parameters defined for boundary faces

real(dp), dimension(:), allocatable :: srds,dns    ! srds = |are|/|dns|, dns = normal distance to cell center from face |dpn*face_normal_unit_vec|
real(dp), dimension(:), allocatable :: srdw,dnw    ! srdw = |are|/|dnw|, dnw = normal distance to cell center from face |dpn*face_normal_unit_vec|
real(dp), dimension(:), allocatable :: srdoc       ! srdoc = |are|/|dpn*face_normal_unit_vec| 
real(dp), dimension(:), allocatable :: foc         ! Interpolation factor for faces at boundaries (known as o-c- faces in structured code)
real(dp), dimension(:), allocatable :: fpro        ! Interpolation factor for faces at processor boundaries

! Mesh topology information - connectivity of cells trough faces
integer, dimension(:), allocatable :: owner     ! Index of the face owner cell
integer, dimension(:), allocatable :: neighbour ! Index of the neighbour cell  - it shares the face with owner
integer, dimension(:), allocatable :: ijl, ijr  ! left and right cell at the block boundary, it is still an inner face in the domain
integer, dimension(:), allocatable :: ijlFace, ijrFace  ! left and right cell face at the block boundary, it is still an inner face in the domain


public 

contains


subroutine triangular_face_area_components_polymesh(px,py,pz,qx,qy,qz,nx,ny,nz)
!
! Ai, i=1,n are n triangular faces enclosing polyhedron P.
! Vertices of Ai are (ai,bi,ci),
!                                _    _   _    _     _  _    _  _
! Unit normal to P on each Ai is ni = ni/|ni|, ni = (bi-ai)x(ci-ai)
! and 
! _    _  _    _    _  _
! p = (bi-ai); q = (ci-ai)
!
! Finally:              _
! Area A of Ai, A = 1/2|ni|
!
! Sources: 
! [1] Dr Robert Nurnberg, Imperial College London, www.ma.ic.ac.uk/~rn/centroid.pdf
! [2] paulbourke.net/geometry/polygonmesh
!
  implicit none
  real(dp), intent(in) :: px,py,pz,qx,qy,qz
  real(dp), intent(inout) :: nx,ny,nz

  ! real(dp), parameter :: half = 0.5_dp

  ! Cross products for triangle surface vectors
  nx = py*qz-pz*qy
  ny = pz*qx-px*qz
  nz = px*qy-py*qx

end subroutine


pure function cell_volume_part_polymesh(ax,ay,az,nx,ny,nz) result(volume_part)
!
! Ai, i=1,n are n triangular faces enclosing polyhedron P.
! Vertices of Ai are (ai,bi,ci),
!                                _    _   _    _
! Unit normal to P on each Ai is ni = ni/|ni|, ni = (bi-ai)x(ci-ai)
! Volume of P is given by:
!                         _  _                     _  _ 
!    V = 1/6 * sum_i=1^n (ai.ni) or sum_i=1^n 1/6*(ai.ni)
!
! Sources: 
! [1] Dr Robert Nurnberg, Imperial College London, www.ma.ic.ac.uk/~rn/centroid.pdf
! [2] paulbourke.net/geometry/polygonmesh

  implicit none
  real(dp), intent(in) :: ax,ay,az
  real(dp), intent(in) :: nx,ny,nz
  real(dp) :: volume_part

  volume_part = 1./6._dp*( ax*nx + ay*ny + az*nz )
end function


pure function centroid_component_part_polymesh(ax,bx,cx,nx,vol) result(cc_part)
!
! Ai, i=1,n are n triangular faces enclosing polyhedron P.
! Vertices of Ai are (ai,bi,ci),
!                                _    _   _    _
! Unit normal to P on each Ai is ni = ni/|ni|, ni = (bi-ai)x(ci-ai)
! Cell center (cc) component of P is given by:
!                                  
!    cc.ed = 1/2V * sum_i=1^n (x.ed)^2*(ni.ed), where ed is unit basis of R^3, d=1,2,3.
!
! Thios function calculates part under summation which is accumulated in outer loop. 
!
! Sources: 
! [1] Dr Robert Nurnberg, Imperial College London, www.ma.ic.ac.uk/~rn/centroid.pdf
! [2] paulbourke.net/geometry/polygonmesh

  implicit none
  real(dp), intent(in) :: ax,bx,cx
  real(dp), intent(in) :: vol,nx
  real(dp) :: cc_part

  cc_part = 1./(2*vol) * 1./24._dp*nx * ( (ax+bx)**2 + (bx+cx)**2 + (cx+ax)**2 )

end function  


  function determinant(a1,a2,a3,b1,b2,b3,q1,q2,q3)
!
! Calculates determinant of 3x3 matrix
!
  implicit none
  real(dp) :: a1,a2,a3,b1,b2,b3,q1,q2,q3
  real(dp) :: determinant

  determinant = (a2*b3-b2*a3)*q1 &
               +(b1*a3-a1*b3)*q2 &
               +(a1*b2-a2*b1)*q3
  end function



subroutine find_intersection_point( &
!                     plane defined by three face corners:
                       x1,y1,z1,&
                       x2,y2,z2, &
                       x3,y3,z3, &
!                      line defined by cell center and neighbour center:
                       x4,y4,z4, &
                       x5,y5,z5, &
!                      intersection point (output):
                       xjp,yjp,zjp &
                       )
!
!***********************************************************************
! Find intersection point (pjp={xjp,yjp,zjp}) of 
! plane (defined by points p1={x1,y1,z1}, p2={x2,y2,z2} and p3={x3,y3,z3}),
! and line (defined by points p4={x4,y4,z4} and p5={x5,y5,z5}).
! The intersection point j' is not the face center j on non-orthogonal meshes.
! There is an "intersection point offset" |jj'| which determines the level
! of nonorthogonality.
!
!
!       |1  1  1  1 |     |1  1  1  0    |
! t = - |x1 x2 x3 x4|  /  |x1 x2 x3 x5-x4|  (mind the minus sign!)
!       |y1 y2 y3 y4| /   |y1 y2 y3 y5-y4|
!       |z1 z2 z3 z4|     |z1 z2 z3 z5-z4|
!
! And intersection point is given by:
! xj = x4 +(x5-x4)*t
! yj = y4 +(y5-y4)*t
! zj = z4 +(z5-z4)*t
!
!
! Nikola Mirkov @2016
!
! example usage: 
! call find_intersection_point( &
!!                            plane defined by face corners:
!                             x(inp),y(inp),z(inp),&
!                             x(inp-idns),y(inp-idns),z(inp-idns), &
!                             x(inp-idtb),y(inp-idtb),z(inp-idtb), &
!!                            line defined by cell center and neighbour center:
!                             xc(inp),yc(inp),zc(inp), &
!                             xc(inp+idew),yc(inp+idew),zc(inp+idew), &
!!                            intersection point:
!                             xj,yj,zj &
!                             )
!***********************************************************************
  implicit none 
!
!***********************************************************************
!

  real(dp), intent(in) :: x1,y1,z1,&
                          x2,y2,z2, &
                          x3,y3,z3, &
                          x4,y4,z4, &
                          x5,y5,z5
  real(dp), intent(inout) :: xjp,yjp,zjp

  real(dp) :: t

 ! Produced by MATLAB symbolic tool.
 t =-(x2*(y3*z4-y4*z3)-x1*(y3*z4-y4*z3)-x3*(y2*z4-y4*z2)+x1*(y2*z4-y4*z2)+x3*(y1*z4-y4*z1)-x2* &
     (y1*z4-y4*z1)+x4*(y2*z3-y3*z2)-x1*(y2*z3-y3*z2)-x4*(y1*z3-y3*z1)+x2*(y1*z3-y3*z1)+x4* &
     (y1*z2-y2*z1)-x3*(y1*z2-y2*z1)) &
    /(x2*(y3*(z5-z4)-(y5-y4)*z3)-x1*(y3*(z5-z4)-(y5-y4)*z3)-x3*(y2*(z5-z4)-(y5-y4)*z2)+x1* &
     (y2*(z5-z4)-(y5-y4)*z2)+x3*(y1*(z5-z4)-(y5-y4)*z1)-x2*(y1*(z5-z4)-(y5-y4)*z1)+(x5-x4)* &
     (y2*z3-y3*z2)-(x5-x4)*(y1*z3-y3*z1)+(x5-x4)*(y1*z2-y2*z1) + tiny)

  xjp = x4 +(x5-x4)*t
  yjp = y4 +(y5-y4)*t
  zjp = z4 +(z5-z4)*t

end subroutine


subroutine read_line_faces_file_polyMesh(faces_file,nn,nod,nmax)
  implicit none
  integer, intent(in) :: faces_file
  integer, intent(in) :: nmax
  integer, intent(out) :: nn
  integer, dimension(nmax), intent(out) :: nod
  integer :: j,m,n
  character(len=15) :: char_string,char_string2

    nn = 0
    nod = 0

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
      nod(1) = m + 1
      nod(2:j-1) = nod(2:j-1) + 1
      nod(nn) = n + 1

end subroutine



subroutine mesh_geometry
!
!  Description:
!    Calculates basic geometrical quantities of numerical mesh
!    defined in this module and needed for FVM computation.
!    Mesh is described in files similar to polyMesh format.
!    Mesh files are 'points', 'faces', 'owner, 'neighbour', 'boundary'
!
!  Date:
!    26/11/2015
!
!  Author:
!    Nikola Mirkov nmirkov@vinca.rs
!
 
  use my_mpi_module

  implicit none

  ! Locals
  integer :: i,j,k,l,ioc,icyc
  integer :: iface
  integer :: inp,inn
  integer :: input_status

  character(len=1) :: ch
  character(len=15) :: char_string,char_string2
  character(len=80) :: line_string
  character(len=10) :: bctype
  character( len = 5) :: nproc_char

  integer, dimension(nomax) :: node  ! It will store global node numbers of cell vertexes
  integer :: nnodes                  ! no. of nodes in face
  integer :: inode                   ! int counter

  integer :: numBoundaries,nfaces,startFace
  integer :: ifaceFriend,startFaceFriend

  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: one_third = 1._dp/3._dp

  real(dp) :: px,py,pz, qx,qy,qz, ax,ay,az, nx,ny,nz, cx,cy,cz
  real(dp) :: xpn,ypn,zpn
  real(dp) :: xjp,yjp,zjp
  real(dp) :: dpn,djn,are
 

!
! > OPEN polyMesh format files: 'points', 'faces', 'owner', 'neighbour'.
!

  ! nproc_char <- myid zapisan levo u vidu stringa.
  call i4_to_s_left ( myid, nproc_char )


  call get_unit( points_file )
  open( unit = points_file,file = 'processor'//trim(nproc_char)//'/constant/polyMesh/points' )
  rewind points_file

  call get_unit( faces_file )
  open( unit = faces_file, file = 'processor'//trim(nproc_char)//'/constant/polyMesh/faces' )
  rewind faces_file

  call get_unit( owner_file )
  open( unit = owner_file, file = 'processor'//trim(nproc_char)//'/constant/polyMesh/owner' )
  rewind owner_file

  call get_unit( neighbour_file )
  open( unit = neighbour_file, file = 'processor'//trim(nproc_char)//'/constant/polyMesh/neighbour' )
  rewind neighbour_file

  call get_unit( boundary_file )
  open( unit = boundary_file, file = 'processor'//trim(nproc_char)//'/constant/polyMesh/boundary' )
  rewind boundary_file

  call get_unit( process_file )
  open( unit = process_file, file = 'processor'//trim(nproc_char)//'/constant/polyMesh/process' )
  rewind process_file


  native_mesh_files = .true.

  ! First characters in OpenFOAM polyMesh files is comment '/*'
  read(points_file,'(a)') ch
  backspace(points_file)
  if(ch == '/')  native_mesh_files = .false.


!
! > Find out number of faces belonging to a specific boundary condition.
! 
  ! Initialize number of boundary faces for each boundary type 
  ninl = 0 
  nout = 0
  nsym = 0
  nwal = 0
  npru = 0
  noc  = 0
  ncyc = 0
  npro = 0

  ! Different types of walls: isothermal, adiabatic, flux.
  nwali = 0
  nwala = 0
  nwalf = 0

  ! Initialize no. of O-C- boundaries
  iOCBoundaries  = 0

  ! Initialize no. of cyclic boundaries
  iCycBoundaries  = 0

  ! Number of rows in the file excluding #comment in header
  call file_row_count ( boundary_file, numBoundaries )

  read(boundary_file,'(a)') line_string
  ! read(boundary_file,*) numBoundaries

  do i=1,numBoundaries
  read(boundary_file,*) bctype,nfaces,startFace

    select case ( adjustl(bctype) )

      case ('inlet')
        if (ninl==0) iInletFacesStart = startFace
        ninl = ninl + nfaces

      case ('outlet')
        if (nout==0) iOutletFacesStart = startFace
        nout = nout + nfaces

      case ('symmetry')
        if (nsym==0) iSymmetryFacesStart = startFace
        nsym = nsym + nfaces

      case ('wall')
        if (nwal==0) iWallFacesStart = startFace
        nwal = nwal + nfaces

      case ('wallIsoth')
        if (nwal==0) iWallFacesStart = startFace
        nwal = nwal + nfaces
        nwali = nwali + nfaces

      case ('wallAdiab')
        if (nwal==0) iWallFacesStart = startFace
        nwal = nwal + nfaces
        nwala = nwala + nfaces

      case ('wallQFlux')
        if (nwal==0) iWallFacesStart = startFace
        nwal = nwal + nfaces
        nwalf = nwalf + nfaces

      case ('prOutlet')
        if (npru==0) iPressOutletFacesStart = startFace
        npru = npru + nfaces

      case ('domain')
        if (noc==0) iOCFacesStart = startFace
        iOCBoundaries = iOCBoundaries + 1
        noc = noc + nfaces

      case ('cyclic')
        if (ncyc==0) iCycFacesStart = startFace
        iCycBoundaries = iCycBoundaries + 1
        ncyc = ncyc + nfaces

      case default
        write(*,*) "Non-existing boundary type in polymesh/boundary file!"
        stop

    end select 

  enddo  

!
! > Read 'process' file with domain connectivity information for MPI
!

  ! Number of rows in the file excluding #comment in header
  call file_row_count ( process_file, numConnections)

  read(process_file,'(a)') line_string
  ! read(process_file,*) numConnections
  
  ! Allocate domain connectivity arrays
  allocate ( neighbProcNo( numConnections ) )
  allocate ( neighbProcOffset( numConnections+1 ) )

  neighbProcOffset(1) = 1

  do i=1,numConnections
    read(process_file,*) neighbProcNo(i),nfaces,startFace
        if (npro==0) iProcFacesStart = startFace
        npro = npro + nfaces
        neighbProcOffset(i+1) = neighbProcOffset(i) + nfaces
  enddo

  ! This denotes position of the last plus one. 
  ! neighbProcOffset( numConnections+1 ) = npro + 1


!
! > Find out numNodes, numFaces, numInnerFaces, etc.
! 

if (native_mesh_files)  then

    ! 'owner' file
    read(owner_file,*) numNodes
    read(owner_file,*) numCells
    read(owner_file,*) numFaces
    read(owner_file,*) numInnerFaces


  else 
  !
  !...if OpenFOAM polyMesh format- the essence of what can be found in polyMeshToCappuccino code.
  ! Code here is tested for OpenFOAM version 4.0, if they don't change polyMesh format this should be OK.
  !       


    ! 'points' file

    ! do i=1,16
    !   read(points_file,*) ch
    ! end do

    do
      read(points_file,*) ch
      if (ch == "(") then
        ! Return two lines
        backspace(points_file)
        backspace(points_file)
        exit
      endif
    end do

    read(points_file,*) numNodes
    read(points_file,*) ch ! reads "("

    
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


    ! 'neighbour' file

    do i=1,17
      read(neighbour_file,*) ch
    end do

    ! do
    !   read(neighbour_file,*) ch
    !   if (ch == "(") then
    !     ! Return two lines
    !     backspace(neighbour_file)
    !     backspace(neighbour_file)
    !     exit
    !   endif
    ! end do

    read(neighbour_file,*) numInnerFaces
    read(neighbour_file,*) ch ! reads "("


    ! 'faces' file

    do i=1,18
      read(faces_file,*) ch
    end do

    ! do
    !   read(faces_file,*) ch
    !   if (ch == "(") then
    !     ! Return two lines
    !     ! backspace(faces_file)
    !     ! backspace(faces_file)
    !     exit
    !   endif
    ! end do

    ! read(faces_file, *) numFaces
    ! read(faces_file,*) ch ! reads "("


  endif

  ! IMPORTANT:
  ! Number of non-zero elements in sparse matrix: nnz
  nnz = 2*numInnerFaces + numCells

  ! IMPORTANT:
  ! Number of boundary faces
  numBoundaryFaces = numFaces - numInnerFaces

  ! IMPORTANT:
  ! Size of arrays storing variables numCells+numBoundaryFaces
  numTotal = numCells + numBoundaryFaces

  ! Size of array for gradients etc. which are numCells plus No. of buffer cells npro
  numPCells = numCells + npro

  ! IMPORTANT:
  ! Where in variable arrays, the boundary face values are stored, 
  ! e.g. for inlet faces it is: (iInletStart+1, iInletStart+Ninl),
  ! for outlet faces it is: (iOutletStart+1, iOutletStart+Nout) etc.
  ! Having all values of a variable, incuding those of volume field (defined in cell centers)
  ! and those of boundary surface field in a single array is helpful.
  iProcStart = numCells

  iInletStart = numCells+Npro

  iOutletStart = numCells+Npro+Ninl

  iSymmetryStart = numCells+Npro+Ninl+Nout

  iWallStart  = numCells+Npro+Ninl+Nout+Nsym

  iPressOutletStart = numCells+Npro+Ninl+Nout+Nsym+Nwal
  
  iOCStart = numCells+Npro+Ninl+Nout+Nsym+Nwal+Npru

  iCycStart = numCells+Npro+Ninl+Nout+Nsym+Nwal+Npru+Noc


!
! > Write report on mesh size into log file
!

  gloNodes = numNodes
  call global_isum( gloNodes )

  gloCells = numCells
  call global_isum( gloCells )

  gloFaces = numFaces
  call global_isum( gloFaces )

  gloIFaces = numInnerFaces
  call global_isum( gloIFaces )
 
  gloBFaces = numBoundaryFaces
  call global_isum( gloBFaces )  

  if ( myid .eq. 0) then

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Mesh data: '



  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nodes, numNodes = ', gloNodes


  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of cells, numCells = ', gloCells


  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of cell-faces, numFaces = ', gloFaces


  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of inner cell-faces, numInnerFaces = ', gloIFaces


  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Boundary information:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of cell-faces on boundary, numBoundaryFaces = ', gloBFaces


  if( ninl.gt.0 ) then 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of inlet faces  = ', ninl
  endif

  if( nout.gt.0 ) then 
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of outlet faces  = ', nout
  endif

  if( nsym.gt.0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of symmetry faces  = ', nsym
  endif

  if( nwal.gt.0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of wall faces  = ', nwal
    if( nwali.gt.0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a,i8)' ) 'Number of isothermal wall faces  = ', nwali
    endif
    if( nwala.gt.0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a,i8)' ) 'Number of adiabatic wall faces  = ', nwala
    endif
    if( nwalf.gt.0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(4x,a,i8)' ) 'Number of flux wall faces  = ', nwalf
    endif
  endif

  if( npru.gt.0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of pressure-outlet faces  = ', npru
  endif
    
  if( noc.gt.0 ) then
    ! Make one correction - divide noc by two because we had read both 'left' and 'right' faces.
    noc = noc/2
    iOCBoundaries = iOCBoundaries/2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of O-C- boundaries  = ', iOCBoundaries
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of O-C- faces  = ', noc
  endif

  if( ncyc.gt.0 ) then
    ! Make one correction - divide ncyc by two because we had read both 'left' and 'right' faces.
    ncyc = ncyc/2
    iCycBoundaries = iCycBoundaries/2
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of cyclic boundaries  = ', iCycBoundaries
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of cyclic faces  = ', ncyc
  endif

  if( npro.gt.0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Number of processor boundary faces  = ', npro
  endif

endif


!
! > Allocate arrays for Mesh description
!

  allocate ( x(numNodes) )
  allocate ( y(numNodes) )
  allocate ( z(numNodes) )

  allocate ( xc(numCells+Npro) )
  allocate ( yc(numCells+Npro) )
  allocate ( zc(numCells+Npro) )

  allocate ( vol(numCells+Npro) )

  allocate ( wallDistance(numCells) )

  allocate ( arx(numFaces) )
  allocate ( ary(numFaces) )
  allocate ( arz(numFaces) )

  allocate ( xf(numFaces) )
  allocate ( yf(numFaces) ) 
  allocate ( zf(numFaces) )

  allocate ( facint(numInnerFaces) ) 

  allocate ( owner(numFaces) )
  allocate ( neighbour(numInnerFaces) )

  allocate ( dns(nsym) )      
  allocate ( srds(nsym) )  

  allocate ( dnw(nwal) )      
  allocate ( srdw(nwal) ) 

  allocate ( ijl(noc+ncyc) ) 
  allocate ( ijr(noc+ncyc) ) 
  allocate ( ijlFace(noc+ncyc) ) 
  allocate ( ijrFace(noc+ncyc) ) 

  allocate ( srdoc(noc+ncyc) ) 
  allocate ( foc(noc+ncyc) )

  allocate ( fpro(npro) )                               

!
! > Allocate and initialize parameters for MPI communication
!
  lenbuf = npro
  allocate ( bufind(lenbuf) )
  allocate ( buffer(lenbuf) )

!
! > Read and process Mesh files, fill owner, neighbour arrays
!

  if (native_mesh_files) then
    
    ! 'points' file
    do i=1,numNodes
        read( points_file, * ) x(i), y(i), z(i) 
    enddo

    ! 'owner' file
    do iface=1,numFaces
      read( owner_file, * ) owner(iface)
    enddo

    ! 'neighbour' file
    do iface=1,numInnerFaces
      read( neighbour_file, * ) neighbour(iface)
    enddo

  else ! OpenFOAM polyMesh format

    ! 'points' file
    do i=1,numNodes
      read(points_file,*) char_string,y(i),char_string2

      ! Char to double conversion:
      read(char_string(2:),*) x(i)
      read(char_string2(1:len_trim(char_string2)-1),*) z(i)  
    end do


    ! 'owner' file
    do i=1,numFaces
      read(owner_file,*) owner(i)
      owner(i) = owner(i) + 1 ! fortran starts from 1
    end do



    ! 'neighbour' file
    do i=1,numInnerFaces
      read(neighbour_file,*) neighbour(i)
      neighbour(i) = neighbour(i) + 1 ! fortran starts from 1
    end do

  endif


  ! Array of buffer cell indexes for MPI exchange
  do i=1,npro
    bufind(i) = owner(iProcFacesStart+i)
  enddo


  !
  ! Rewind boundary file for one more sweep - to read O-C- and cyclic boundaries
  !
  if(noc.gt.0 .or. ncyc.gt.0) then

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Defining O-C- and cyclic boundary face pairs...'

    rewind( boundary_file )

    read(boundary_file,'(a)') line_string
    ! read(boundary_file,*) numBoundaries

    ioc = 1
    icyc = 1

    bc_loop: do

      read(boundary_file,*,iostat = input_status) bctype,nfaces,startFace 

      if(input_status /= 0) exit


      if ( adjustl(bctype) == 'domain' ) then ! we have what we call a o-c- boundary

        ! We expect to read line with info where to find corresponding faces (bctype='domain' and nfaces is the same of course.)
        read(boundary_file,*)  bctype,nfaces,startFaceFriend

        do k=1,nfaces
          iface = startFace + k
          ifaceFriend = startFaceFriend + k

          ijl(ioc) = owner(iface)
          ijr(ioc) = owner(ifaceFriend)

          ijlFace(ioc) = iface
          ijrFace(ioc) = ifaceFriend

          ioc = ioc + 1

        enddo

      elseif ( adjustl(bctype) == 'cyclic' ) then ! we have cyclic boundary

        ! We expect to read line with info where to find corresponding faces (bctype='cyclic' and nfaces is the same of course.)
        read(boundary_file,*)  bctype,nfaces,startFaceFriend

        do k=1,nfaces
          iface = startFace + k
          ifaceFriend = startFaceFriend + k

          ijl(noc+icyc) = owner(iface)
          ijr(noc+icyc) = owner(ifaceFriend)

          ijlFace(noc+icyc) = iface
          ijrFace(noc+icyc) = ifaceFriend

          icyc = icyc + 1

        enddo

      endif
    
    enddo bc_loop

  endif

  !
  ! > Cell volumes, cell face centers
  !

  do iface=1,numFaces

    inp = owner(iface)

    node = 0

    ! Read line in 'faces' file
    if (native_mesh_files) then
      read( faces_file, * ) nnodes, (node(k), k=1,nnodes)
    else ! OpenFOAM polyMesh
      call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)
    endif

    ax = 0.0_dp
    ay = 0.0_dp
    az = 0.0_dp

    do i=1,nnodes-2 
      ! Vectors to vertices
      ! 2-1
      px = x(node(i+1))-x(node(1))
      py = y(node(i+1))-y(node(1))
      pz = z(node(i+1))-z(node(1))
      ! 3-1
      qx = x(node(i+2))-x(node(1))
      qy = y(node(i+2))-y(node(1))
      qz = z(node(i+2))-z(node(1))

      call triangular_face_area_components_polymesh( px,py,pz, qx,qy,qz, nx,ny,nz )
 
      !
      ! > Cell-face area vector components (Area vector lies in direction of face normal)
      ! 

      arx(iface) = arx(iface) + half*nx
      ary(iface) = ary(iface) + half*ny
      arz(iface) = arz(iface) + half*nz
 
      !
      ! > Cell-face centroid components - accumulation stage
      !
      
      cx = one_third*( x(node(i+2)) + x(node(i+1)) + x(node(1)) )
      cy = one_third*( y(node(i+2)) + y(node(i+1)) + y(node(1)) )
      cz = one_third*( z(node(i+2)) + z(node(i+1)) + z(node(1)) )

      xf(iface) = xf(iface) + abs(nx*cx)
      yf(iface) = yf(iface) + abs(ny*cy)
      zf(iface) = zf(iface) + abs(nz*cz)

      ax = ax + abs(nx)
      ay = ay + abs(ny)
      az = az + abs(nz)

      !
      ! > Compute cell volumes 
      !

      vol(inp) = vol(inp) + cell_volume_part_polymesh( cx, cy, cz, nx, ny, nz ) 
      if ( iface <= numInnerFaces ) then 
        inn = neighbour(iface)
        vol(inn) = vol(inn) + cell_volume_part_polymesh( cx, cy, cz, -nx, -ny, -nz )
      endif

    enddo

    ! ! > Cell-face centroid components - final
    !  if(iface.le.numInnerFaces) then
    !     xf(iface) = xf(iface) / (ax+tiny)
    !     yf(iface) = yf(iface) / (ay+tiny)
    !     zf(iface) = zf(iface) / (az+tiny)
    ! else   
        ! > Because I could have not resolved the problem, these line are inserted 
        !   where face centroid is calculated by arithmetic average.  
        xf(iface) = 0.0_dp
        yf(iface) = 0.0_dp
        zf(iface) = 0.0_dp

        do inode=1,nnodes
          xf(iface) = xf(iface)+x(node(inode))
          yf(iface) = yf(iface)+y(node(inode))
          zf(iface) = zf(iface)+z(node(inode))
        end do

        xf(iface) = xf(iface) / dble(nnodes)
        yf(iface) = yf(iface) / dble(nnodes)
        zf(iface) = zf(iface) / dble(nnodes)
    ! endif

  enddo


  ! Rewind 'faces' file for one more sweep
  rewind( faces_file )

  if (.not.native_mesh_files) then

    do i=1,18
      read(faces_file,*) ch
    end do

    ! do
    !   read(faces_file,*) ch
    !   if (ch == "(") then
    !     exit
    !   endif
    ! end do  

  endif

  !
  ! > Cell centers
  !

  do iface=1,numFaces

    inp = owner(iface)

    node = 0

    ! Read line in 'faces' file
    if (native_mesh_files) then
      read( faces_file, * ) nnodes, (node(k), k=1,nnodes)
    else ! OpenFOAM polyMesh
      call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)
    endif

    do i=1,nnodes-2 
      ! Vectors to vertices
      ! 2-1
      px = x(node(i+1))-x(node(1))
      py = y(node(i+1))-y(node(1))
      pz = z(node(i+1))-z(node(1))
      ! 3-1
      qx = x(node(i+2))-x(node(1))
      qy = y(node(i+2))-y(node(1))
      qz = z(node(i+2))-z(node(1))
   
      call triangular_face_area_components_polymesh( px,py,pz, qx,qy,qz, nx,ny,nz )

      xc(inp) = xc(inp) + centroid_component_part_polymesh( x(node(1)),x(node(i+1)),x(node(i+2)), nx, vol(inp) )
      yc(inp) = yc(inp) + centroid_component_part_polymesh( y(node(1)),y(node(i+1)),y(node(i+2)), ny, vol(inp) )
      zc(inp) = zc(inp) + centroid_component_part_polymesh( z(node(1)),z(node(i+1)),z(node(i+2)), nz, vol(inp) )

      if ( iface <= numInnerFaces ) then
      inn = neighbour(iface)
      xc(inn) = xc(inn) + centroid_component_part_polymesh( x(node(1)),x(node(i+2)),x(node(i+1)), -nx, vol(inn) )
      yc(inn) = yc(inn) + centroid_component_part_polymesh( y(node(1)),y(node(i+2)),y(node(i+1)), -ny, vol(inn) )
      zc(inn) = zc(inn) + centroid_component_part_polymesh( z(node(1)),z(node(i+2)),z(node(i+1)), -nz, vol(inn) )
      endif

    enddo

  enddo


  ! We need some geometry in the buffer cells, we fill these by exchanging info with other processes
  call exchange( xc )
  call exchange( yc )
  call exchange( zc )
  call exchange( Vol )


  ! Rewind 'faces' file for one more sweep
  rewind( faces_file )

  if (.not.native_mesh_files) then

    do i=1,18
      read(faces_file,*) ch
    end do

    ! do
    !   read(faces_file,*) ch
    !   if (ch == "(") then
    !     exit
    !   endif
    ! end do 

  endif

  !
  ! > Interpolation factor > inner faces
  !

  do iface=1,numInnerFaces

    inp = owner(iface)
    inn = neighbour(iface)

    node = 0

    ! Read line in 'faces' file
    if (native_mesh_files) then
      read( faces_file, * ) nnodes, (node(k), k=1,nnodes)
    else ! OpenFOAM polyMesh
      call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)
    endif

    xpn = xc(inn)-xc(inp)
    ypn = yc(inn)-yc(inp)
    zpn = zc(inn)-zc(inp)

    dpn = sqrt( xpn**2 + ypn**2 + zpn**2 )

    ! > > Intersection point j' of line connecting centers with cell face, we are taking only three points assuming that other are co-planar
    call find_intersection_point( &
                                 ! plane defined by three face vertices:
                                 x(node(1)),y(node(1)),z(node(1)),&
                                 x(node(2)),y(node(2)),z(node(2)), &
                                 x(node(3)),y(node(3)),z(node(3)), &
                                 ! line defined by cell center and neighbour center:
                                 xc(inp),yc(inp),zc(inp), &
                                 xc(inn),yc(inn),zc(inn), &
                                 ! intersection point (output):
                                 xjp,yjp,zjp &
                                )
    xpn = xjp - xc(inp)
    ypn = yjp - yc(inp)
    zpn = zjp - zc(inp)

    djn = sqrt( xpn**2 + ypn**2 + zpn**2 )

    ! Interpolation factor |P Pj'|/|P Pj| where P is cell center, Pj neighbour cell center and j' intersection point.
    facint(iface) = djn/dpn
    
  enddo


  !
  ! > Interpolation factor > faces on process boundaries
  !

  ! Rewind 'faces' file for one more sweep
  rewind( faces_file )
  if (.not.native_mesh_files) then

    do i=1,18
      read(faces_file,*) ch
    end do

    ! do
    !   read(faces_file,*) ch
    !   if (ch == "(") then
    !     exit
    !   endif
    ! end do 

  endif
  
  ! Fast forward to place where processor faces start
  do i=1,iProcFacesStart
    read(faces_file,*) ch
  enddo

  do i=1,npro
    iface = iProcFacesStart + i
    inp = owner( iface )
    inn = iProcStart + i

    node = 0

    ! Read line in 'faces' file
    if (native_mesh_files) then
      read( faces_file, * ) nnodes, (node(k), k=1,nnodes)
    else ! OpenFOAM polyMesh
      call read_line_faces_file_polyMesh(faces_file,nnodes,node,nomax)
    endif

    xpn = xc(inn)-xc(inp)
    ypn = yc(inn)-yc(inp)
    zpn = zc(inn)-zc(inp)

    dpn = sqrt( xpn**2 + ypn**2 + zpn**2 )

    ! > > Intersection point j' of line connecting centers with cell face, we are taking only three points assuming that other are co-planar
    call find_intersection_point( &
                                 ! plane defined by three face vertices:
                                 x(node(1)),y(node(1)),z(node(1)),&
                                 x(node(2)),y(node(2)),z(node(2)), &
                                 x(node(3)),y(node(3)),z(node(3)), &
                                 ! line defined by cell center and neighbour center:
                                 xc(inp),yc(inp),zc(inp), &
                                 xc(inn),yc(inn),zc(inn), &
                                 ! intersection point (output):
                                 xjp,yjp,zjp &
                                )
    xpn = xjp - xc(inp)
    ypn = yjp - yc(inp)
    zpn = zjp - zc(inp)

    djn = sqrt( xpn**2 + ypn**2 + zpn**2 )

    ! Interpolation factor |P Pj'|/|P Pj| where P is cell center, Pj neighbour cell center and j' intersection point.
    fpro(i) = djn/dpn

  enddo


  !
  ! > Interpolation factor > O-C boundaries
  !
  if(noc.gt.0) then 
    do ioc=1,noc
          iface = ijlFace(ioc)
          ifaceFriend = ijrFace(ioc)

          inp = ijl(ioc)
          inn = ijr(ioc)


          xpn = xc(inn) - xf( ifaceFriend )
          ypn = yc(inn) - yf( ifaceFriend )
          zpn = zc(inn) - zf( ifaceFriend )

          dpn = sqrt( xpn**2 + ypn**2 + zpn**2 )


          xpn = xf(iface) - xc(inp)
          ypn = yf(iface) - yc(inp)
          zpn = zf(iface) - zc(inp)

          djn = sqrt( xpn**2 + ypn**2 + zpn**2 )

          ! Interpolation factor
          foc(ioc) = djn/(djn+dpn)

          ! Face area 
          are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

          ! Cell face area divided by distance to the cell center
          srdoc(ioc) = are/(djn+dpn)

    enddo
  endif

  !
  ! > Interpolation factor > cyclic boundaries
  !
  if(ncyc.gt.0) then 
    do icyc=1,ncyc
          iface = ijlFace(noc+icyc)
          ifaceFriend = ijrFace(noc+icyc)

          inp = ijl(noc+icyc)
          inn = ijr(noc+icyc)


          xpn = xc(inn) - xf( ifaceFriend )
          ypn = yc(inn) - yf( ifaceFriend )
          zpn = zc(inn) - zf( ifaceFriend )

          dpn = sqrt( xpn**2 + ypn**2 + zpn**2 )


          xpn = xf(iface) - xc(inp)
          ypn = yf(iface) - yc(inp)
          zpn = zf(iface) - zc(inp)

          djn = sqrt( xpn**2 + ypn**2 + zpn**2 )

          ! ! Interpolation factor
          ! foc(noc+icyc) = djn/(djn+dpn)

          ! U sustini uvek je periodic pravac sa uniformnom mrezom pa je:
          foc(noc+icyc) = 0.5_dp

          ! Face area 
          are = sqrt(arx(iface)**2+ary(iface)**2+arz(iface)**2)

          ! Cell face area divided by distance to the cell center
          srdoc(noc+icyc) = are/(djn+dpn)

    enddo
  endif
  
  ! Important!
  ! I would like to make it simple for the rest of the code so the massive noc+ncyc
  ! is encapsulated in  single number.
  noc = noc + ncyc


!
! > Report on geometrical quantities
!

  if (myid .eq. 0 ) then

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cell data: '

  call r8vec_print_some ( numCells, vol, 1, 10, &
      '  First 10 elements of cell volumes array:' )

  call r8vec_print_some ( numCells, xc, 1, 10, &
      '  First 10 elements of cell x-centers array:' )

  call r8vec_print_some ( numCells, yc, 1, 10, &
      '  First 10 elements of cell y-centers array:' )

  call r8vec_print_some ( numCells, zc, 1, 10, &
      '  First 10 elements of cell z-centers array:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Face data: '

  call i4vec_print2 ( 10, owner, neighbour, '  First 10 lines of owner and neighbour arrays:' )

  call r8vec_print_some ( numFaces, arx, 1, 10, &
      '  First 10 elements of Arx array:' )

  call r8vec_print_some ( numFaces, ary, 1, 10, &
      '  First 10 elements of Ary array:' )

  call r8vec_print_some ( numFaces, arz, 1, 10, &
      '  First 10 elements of Arz array:' )

    call r8vec_print_some ( numFaces, xf, 1, 10, &
      '  First 10 elements of xf array:' )

  call r8vec_print_some ( numFaces, yf, 1, 10, &
      '  First 10 elements of yf array:' )

  call r8vec_print_some ( numFaces, zf, 1, 10, &
      '  First 10 elements of zf array:' )

  call r8vec_print_some ( numInnerFaces, facint, 1, 10, &
      '  First 10 elements of interpolation factor (facint) array:' )

  endif


!
!  > CLOSE polyMesh format file: 'points', 'faces', 'owner', 'neighbour', 'boundary', 'process'.
!
  close ( points_file )
  close ( faces_file )
  close ( owner_file )
  close ( neighbour_file)
  close ( boundary_file)
  close ( process_file)
!+-----------------------------------------------------------------------------+


end subroutine mesh_geometry


end module
