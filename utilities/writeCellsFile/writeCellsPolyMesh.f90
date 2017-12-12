!*******************************************************************************
program writeCells
!
!  Description:
!    Creates Cells file in polyMesh directory based on node numbering in
!    blockMesh code.
!
!  Date:
!    11/10/2017
!
!  Author:
!    Nikola Mirkov nmirkov@vinca.rs
!
!*******************************************************************************
!
  use cells_module

  implicit none

  integer, parameter :: hexType = 5
  integer :: cells_file, input_file
  integer :: i,j,k
  integer :: icell,jcell,kcell,icurrent,floor,hallway


  call get_unit( cells_file )
  open( unit = cells_file, file='polyMesh/cells' )
  rewind cells_file

  call get_unit( input_file )
  open( unit = input_file, file='polyMesh/cells_input' )
  rewind input_file

  ! Read number of cells in each direction
  read(input_file,*) icell,jcell,kcell

  do k=1,kcell
    do j=1,jcell
      do i=1,icell
        floor = (icell+1)*(jcell+1)
        hallway = (icell+1)
        icurrent = floor*(k-1) + hallway*(j-1) + i
        write(cells_file,'(9(i0,1x))') hexType, icurrent, icurrent+1, icurrent+icell+2, icurrent+icell+1, &
      &                    icurrent+floor,icurrent+floor+1,icurrent+floor+icell+2,icurrent+floor+icell+1
      enddo
    enddo
  enddo


 close(cells_file)
 close(input_file)

end program