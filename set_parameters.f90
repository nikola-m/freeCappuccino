!***********************************************************************
!
subroutine set_parameters
!
!***********************************************************************
!
!  Mesh parameters
! 
!  Discussion:
!    owner array's lehgth is equal to numFacesTotal, because it contains
!    both inner face owner cells, and boundary face owner cells.
!    
!    neighbour array length is equal to numInnerFaces
!
!    number of non-zero elements in sparse matrix nnz = numCells+2*numInnerFaces
!
!    numTotal will be the length of arrays holding field variables, because
!    it is convinient to store volume field values (defined in cell centers)
!    and those of boundary surface field values in a single array.
!
!***********************************************************************
!
  use parameters
  use geometry, only: numTotal,numCells,numInnerFaces,Ninl,Nout,Nsym,Nwal,Npru,Noc,&
    iInletStart,iOutletStart,iSymmetryStart,iWallStart,iPressOutletStart,iOCStart,&
    iInletFacesStart,iOutletFacesStart,iSymmetryFacesStart,iWallFacesStart,iPressOutletFacesStart,iOCFacesStart

  implicit none
!
!***********************************************************************
!
  numTotal = numCells+Ninl+Nout+Nsym+Nwal+Npru+Noc

  ! numFacesTotal = numInnerFaces+Ninl+Nout+Nsym+Nwal+Npru ! maybe good for checking only..?

  ! nnz = numCells + 2*numInnerFaces

  ! Where in variable arrays, the boundary face values are stored, 
  ! e.g. for inlet faces it is: (iInletStart+1, iInletStart+Ninl),
  ! for outlet faces it is: (iOutletStart+1, iOutletStart+Nout) etc.
  ! Having all values of a variable, incuding those of volume field (defined in cell centers)
  ! and those of boundary surface field in a single array is helpful.
  iInletStart = numCells
  iOutletStart = numCells+Ninl
  iSymmetryStart = numCells+Ninl+Nout
  iWallStart = numCells+Ninl+Nout+Nsym
  iPressOutletStart = numCells+Ninl+Nout+Nsym+Nwal
  iOCStart = numCells+Ninl+Nout+Nsym+Nwal+Npru

  ! The owner array is longer than the neighbour array because it contains the owner cells of
  ! the boundary faces, belonging to different boundary conditions categorised by order.
  ! Therefore, we can find owner cells of inlet faces in positions owner( iInletFacesStart+1 : iInletFacesStart+Ninl ),
  ! wall face owners in positions owner( iWallFacesStart+1 : iWallFacesStart+Ninl ), 
  ! etc.
  iInletFacesStart = numInnerFaces
  iOutletFacesStart = numInnerFaces+Ninl
  iSymmetryFacesStart = numInnerFaces+Ninl+Nout
  iWallFacesStart = numInnerFaces+Ninl+Nout+Nsym
  iPressOutletFacesStart = numInnerFaces+Ninl+Nout+Nsym+Nwal
  iOCFacesStart = numInnerFaces+Ninl+Nout+Nsym+Nwal+Npru

end subroutine
