subroutine grad_lsq(fi,dFidxi,istage,dmat)
!
!***********************************************************************
!
!      Purpose:
!      Calculates cell-centered gradients using UNWEIGHTED Least-Squares approach.
!
!      Description:
!      Approach taken from PhD thesis of Bojan Niceno, TU Delft, 2000.
!
!      Arguments:
!
!      FI - field variable which gradient we look for.
!      DFIDX,DFIDY,DFIDZ - cell centered gradient - a three component gradient vector.
!      ISTAGE - integer. If ISTAGE=1 calculates and stores only geometrical
!      parameters - a system matrix for least square problem at every cell. 
!      Usually it is called with ISTAGE=1 at the beggining of simulation.
!      If 2 it doesn't calculate system matrix, just RHS and solves system.
!      Dmat - LSQ matrix with geometry data
!      XC,YC,ZC - coordinates of cell centers    
!
!      Example call:
!      CALL GRADFI_LSQ(U,DUDX,DUDY,DUDZ,2,D)
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use geometry

  implicit none

  integer, intent(in) :: istage
  real(dp),dimension(numTotal), intent(in)   :: fi
  real(dp),dimension(3,numCells), intent(inout) :: dFidxi
  real(dp),dimension(6,numCells), intent(inout) :: dmat

  !
  ! Locals
  !
  integer :: i,ijp,ijn,inp

  real(dp) :: d11,d12,d13,d21,d22,d23,d31,d32,d33
  real(dp) :: tmp

  real(dp), dimension(numCells) :: b1,b2,b3 
!
!***********************************************************************
!

  ! Initialize dmat matrix:
  dmat(:,:) = 0.0d0

  if(istage.eq.1) then
  ! Coefficient matrix - should be calculated only once 

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
        ijp = owner(i)
        ijn = neighbour(i)

        Dmat(1,ijp) = Dmat(1,ijp) + (xc(ijn)-xc(ijp))**2
        Dmat(1,ijn) = Dmat(1,ijn) + (xc(ijp)-xc(ijn))**2 

        Dmat(4,ijp) = Dmat(4,ijp) + (yc(ijn)-yc(ijp))**2
        Dmat(4,ijn) = Dmat(4,ijn) + (yc(ijp)-yc(ijn))**2

        Dmat(6,ijp) = Dmat(6,ijp) + (zc(ijn)-zc(ijp))**2
        Dmat(6,ijn) = Dmat(6,ijn) + (zc(ijp)-zc(ijn))**2  

        Dmat(2,ijp) = Dmat(2,ijp) + (xc(ijn)-xc(ijp))*(yc(ijn)-yc(ijp)) 
        Dmat(2,ijn) = Dmat(2,ijn) + (xc(ijp)-xc(ijn))*(yc(ijp)-yc(ijn))

        Dmat(3,ijp) = Dmat(3,ijp) + (xc(ijn)-xc(ijp))*(zc(ijn)-zc(ijp)) 
        Dmat(3,ijn) = Dmat(3,ijn) + (xc(ijp)-xc(ijn))*(zc(ijp)-zc(ijn))
   
        Dmat(5,ijp) = Dmat(5,ijp) + (yc(ijn)-yc(ijp))*(zc(ijn)-zc(ijp)) 
        Dmat(5,ijn) = Dmat(5,ijn) + (yc(ijp)-yc(ijn))*(zc(ijp)-zc(ijn))                                                                 
  enddo     
 
  
  ! Faces along O-C grid cuts
  do i=1,noc
    ijp = ijl(i)
    ijn = ijr(i)

        Dmat(1,ijp) = Dmat(1,ijp) + (xc(ijn)-xc(ijp))**2
        Dmat(1,ijn) = Dmat(1,ijn) + (xc(ijp)-xc(ijn))**2 

        Dmat(4,ijp) = Dmat(4,ijp) + (yc(ijn)-yc(ijp))**2
        Dmat(4,ijn) = Dmat(4,ijn) + (yc(ijp)-yc(ijn))**2

        Dmat(6,ijp) = Dmat(6,ijp) + (zc(ijn)-zc(ijp))**2
        Dmat(6,ijn) = Dmat(6,ijn) + (zc(ijp)-zc(ijn))**2  

        Dmat(2,ijp) = Dmat(2,ijp) + (xc(ijn)-xc(ijp))*(yc(ijn)-yc(ijp)) 
        Dmat(2,ijn) = Dmat(2,ijn) + (xc(ijp)-xc(ijn))*(yc(ijp)-yc(ijn))

        Dmat(3,ijp) = Dmat(3,ijp) + (xc(ijn)-xc(ijp))*(zc(ijn)-zc(ijp)) 
        Dmat(3,ijn) = Dmat(3,ijn) + (xc(ijp)-xc(ijn))*(zc(ijp)-zc(ijn))
   
        Dmat(5,ijp) = Dmat(5,ijp) + (yc(ijn)-yc(ijp))*(zc(ijn)-zc(ijp)) 
        Dmat(5,ijn) = Dmat(5,ijn) + (yc(ijp)-yc(ijn))*(zc(ijp)-zc(ijn))  

  end do

  
  ! Boundary faces:

  ! Inlet faces
  do i=1,ninl
    ijp = owner(iInletFacesStart+i)
        Dmat(1,ijp) = Dmat(1,ijp) + (xfi(i)-xc(ijp))**2
        Dmat(4,ijp) = Dmat(4,ijp) + (yfi(i)-yc(ijp))**2
        Dmat(6,ijp) = Dmat(6,ijp) + (zfi(i)-zc(ijp))**2
        Dmat(2,ijp) = Dmat(2,ijp) + (xfi(i)-xc(ijp))*(yfi(i)-yc(ijp)) 
        Dmat(3,ijp) = Dmat(3,ijp) + (xfi(i)-xc(ijp))*(zfi(i)-zc(ijp)) 
        Dmat(5,ijp) = Dmat(5,ijp) + (yfi(i)-yc(ijp))*(zfi(i)-zc(ijp)) 
  end do

  ! Outlet faces
  do i=1,nout
    ijp = owner(iOutletFacesStart+i)
        Dmat(1,ijp) = Dmat(1,ijp) + (xfo(i)-xc(ijp))**2
        Dmat(4,ijp) = Dmat(4,ijp) + (yfo(i)-yc(ijp))**2
        Dmat(6,ijp) = Dmat(6,ijp) + (zfo(i)-zc(ijp))**2
        Dmat(2,ijp) = Dmat(2,ijp) + (xfo(i)-xc(ijp))*(yfo(i)-yc(ijp)) 
        Dmat(3,ijp) = Dmat(3,ijp) + (xfo(i)-xc(ijp))*(zfo(i)-zc(ijp)) 
        Dmat(5,ijp) = Dmat(5,ijp) + (yfo(i)-yc(ijp))*(zfo(i)-zc(ijp)) 
  end do

  ! Symmetry faces
  do i=1,nsym
    ijp = owner(iSymmetryFacesStart+i)
        Dmat(1,ijp) = Dmat(1,ijp) + (xfs(i)-xc(ijp))**2
        Dmat(4,ijp) = Dmat(4,ijp) + (yfs(i)-yc(ijp))**2
        Dmat(6,ijp) = Dmat(6,ijp) + (zfs(i)-zc(ijp))**2
        Dmat(2,ijp) = Dmat(2,ijp) + (xfs(i)-xc(ijp))*(yfs(i)-yc(ijp)) 
        Dmat(3,ijp) = Dmat(3,ijp) + (xfs(i)-xc(ijp))*(zfs(i)-zc(ijp))    
        Dmat(5,ijp) = Dmat(5,ijp) + (yfs(i)-yc(ijp))*(zfs(i)-zc(ijp)) 
  end do

  ! Wall faces
  do i=1,nwal
    ijp = owner(iWallFacesStart+i)
        Dmat(1,ijp) = Dmat(1,ijp) + (xfw(i)-xc(ijp))**2
        Dmat(4,ijp) = Dmat(4,ijp) + (yfw(i)-yc(ijp))**2
        Dmat(6,ijp) = Dmat(6,ijp) + (zfw(i)-zc(ijp))**2
        Dmat(2,ijp) = Dmat(2,ijp) + (xfw(i)-xc(ijp))*(yfw(i)-yc(ijp)) 
        Dmat(3,ijp) = Dmat(3,ijp) + (xfw(i)-xc(ijp))*(zfw(i)-zc(ijp))   
        Dmat(5,ijp) = Dmat(5,ijp) + (yfw(i)-yc(ijp))*(zfw(i)-zc(ijp)) 
  end do

  ! Pressure outlet faces
  do i=1,npru
    ijp = owner(iPressOutletFacesStart+i)
        Dmat(1,ijp) = Dmat(1,ijp) + (xfpr(i)-xc(ijp))**2
        Dmat(4,ijp) = Dmat(4,ijp) + (yfpr(i)-yc(ijp))**2
        Dmat(6,ijp) = Dmat(6,ijp) + (zfpr(i)-zc(ijp))**2
        Dmat(2,ijp) = Dmat(2,ijp) + (xfpr(i)-xc(ijp))*(yfpr(i)-yc(ijp)) 
        Dmat(3,ijp) = Dmat(3,ijp) + (xfpr(i)-xc(ijp))*(zfpr(i)-zc(ijp))   
        Dmat(5,ijp) = Dmat(5,ijp) + (yfpr(i)-yc(ijp))*(zfpr(i)-zc(ijp))    
  end do


  elseif(istage.eq.2) then

  ! Initialize rhs vector
  b1(:) = 0.0d0
  b2(:) = 0.0d0
  b3(:) = 0.0d0

  ! Inner faces:

  do i=1,numInnerFaces                                                       
        ijp = owner(i)
        ijn = neighbour(i)

        b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xc(ijn)-xc(ijp)) 
        b1(ijn) = b1(ijn) + (Fi(ijp)-Fi(ijn))*(xc(ijp)-xc(ijn))

        b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yc(ijn)-yc(ijp)) 
        b2(ijn) = b2(ijn) + (Fi(ijp)-Fi(ijn))*(yc(ijp)-yc(ijn))  

        b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zc(ijn)-zc(ijp)) 
        b3(ijn) = b3(ijn) + (Fi(ijp)-Fi(ijn))*(zc(ijp)-zc(ijn))
                                                                                                                         
  enddo     


  
  ! Faces along O-C grid cuts
  do i=1,noc
    ijp = ijl(i)
    ijn = ijr(i)

        b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xc(ijn)-xc(ijp)) 
        b1(ijn) = b1(ijn) + (Fi(ijp)-Fi(ijn))*(xc(ijp)-xc(ijn))

        b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yc(ijn)-yc(ijp)) 
        b2(ijn) = b2(ijn) + (Fi(ijp)-Fi(ijn))*(yc(ijp)-yc(ijn))  

        b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zc(ijn)-zc(ijp)) 
        b3(ijn) = b3(ijn) + (Fi(ijp)-Fi(ijn))*(zc(ijp)-zc(ijn))
  
  end do

  ! Boundary faces:

  ! Inlet faces
  do i=1,ninl
    ijp = owner(iInletFacesStart+i)
    ijn = iInletStart+i

        b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xfi(i)-xc(ijp)) 
        b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yfi(i)-yc(ijp))  
        b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zfi(i)-zc(ijp)) 
  end do

  ! Outlet faces
  do i=1,nout
    ijp = owner(iOutletFacesStart+i)
    ijn = iOutletStart+i

        b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xfo(i)-xc(ijp)) 
        b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yfo(i)-yc(ijp))  
        b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zfo(i)-zc(ijp)) 

  end do

  ! Symmetry faces
  do i=1,nsym
    ijp = owner(iSymmetryFacesStart+i)
    ijn = iSymmetryStart+i

        b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xfs(i)-xc(ijp)) 
        b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yfs(i)-yc(ijp))  
        b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zfs(i)-zc(ijp)) 
  end do

  ! Wall faces
  do i=1,nwal
    ijp = owner(iWallFacesStart+i)
    ijn = iWallStart+i

        b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xfw(i)-xc(ijp)) 
        b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yfw(i)-yc(ijp))  
        b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zfw(i)-zc(ijp)) 
        
  end do

  ! Pressure outlet faces
  do i=1,npru
    ijp = owner(iPressOutletFacesStart+i)
    ijn = iPressOutletStart+i

        b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xfpr(i)-xc(ijp)) 
        b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yfpr(i)-yc(ijp))  
        b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zfpr(i)-zc(ijp))     
  end do

  ! Calculate gradient

  ! Cell loop
  do inp=1,numCells 

        ! Copy from Coefficient matrix 
        d11 = Dmat(1,inp)
        d12 = Dmat(2,inp)
        d13 = Dmat(3,inp)

        d22 = Dmat(4,inp)
        d23 = Dmat(5,inp)
        d33 = Dmat(6,inp)

        ! Symmetric part
        d21 = d12
        d31 = d13
        d32 = d23

        ! Solve system 
        
        tmp = 1./(d11*d22*d33 - d11*d23*d32 - d12*d21*d33 + d12*d23*d31 + d13*d21*d32 - d13*d22*d31)

        dFidxi(1,inp) = ( (b1(inp)*(d22*d33 - d23*d32)) - (b2(inp)*(d21*d33 - d23*d31)) + (b3(inp)*(d21*d32 - d22*d31)) ) * tmp
        dFidxi(2,inp) = ( (b2(inp)*(d11*d33 - d13*d31)) - (b1(inp)*(d12*d33 - d13*d32)) - (b3(inp)*(d11*d32 - d12*d31)) ) * tmp
        dFidxi(3,inp) = ( (b1(inp)*(d12*d23 - d13*d22)) - (b2(inp)*(d11*d23 - d13*d21)) + (b3(inp)*(d11*d22 - d12*d21)) ) * tmp

  enddo
   

  endif 
  
end subroutine

