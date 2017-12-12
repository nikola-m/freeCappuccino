subroutine grad_lsq_dm(fi,dFidxi,istage,dmat)
!
!***********************************************************************
!
!  Purpose:
!  Calculates cell-centered gradients using WEIGHTED Least-Squares approach.
!
!  Description:
!  Approach taken from a paper:
!  Dimitry Mavriplis, "Revisiting the Least Squares procedure for Gradient Reconstruction on Unstructured Meshes." NASA/CR-2003-212683.
!  Weights based on inverse cell-centers distance is added to improve conditioning of the system matrix for skewed meshes.
!  Reduces storage requirements compared to QR subroutine.
!  System matrix is symmetric and can be solved efficiently using Cholesky decomposition or by matrix inversion. 
!
!  Arguments:
!
!  FI - field variable which gradient we look for.
!  DFIDXi- cell centered gradient - a three component gradient vector.
!  ISTAGE - integer. If ISTAGE=1 calculates and stores only geometrical
!  parameters - a system matrix for least square problem at every cell. 
!  Usually it is called with ISTAGE=1 at the beggining of simulation.
!  If 2 it doesn't calculate system matrix, just RHS and solves system.
!  Dmat - LSQ matrix with geometry data
!
!  Example call:
!  CALL GRADFI_LSQ_DM(U,DUDX,DUDY,DUDZ,2,D)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry

  implicit none

  integer, intent(in) :: istage
  real(dp),dimension(numTotal), intent(in)   :: fi
  real(dp),dimension(3,numPCells), intent(inout) :: dFidxi
  real(dp),dimension(9,numCells), intent(inout) :: dmat

  !
  ! Locals
  !
  integer :: i,ijp,ijn,inp,iface

  real(dp) :: w
  real(dp) :: Dx,Dy,Dz
  real(dp) :: d11,d12,d13,d21,d22,d23,d31,d32,d33
  real(dp) :: tmp

  real(dp), dimension(numCells) :: b1,b2,b3 
!
!***********************************************************************
!

  if(istage.eq.1) then
  ! Coefficient matrix - should be calculated only once 

  ! Initialize dmat matrix:
  dmat = 0.0d0

  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
        ijp = owner(i)
        ijn = neighbour(i)

        w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

        Dx = xc(ijn)-xc(ijp)
        Dy = yc(ijn)-yc(ijp)
        Dz = zc(ijp)-zc(ijn)

        Dmat(1,ijp) = Dmat(1,ijp) + w*Dx*Dx
        Dmat(1,ijn) = Dmat(1,ijn) + w*Dx*Dx 

        Dmat(4,ijp) = Dmat(4,ijp) + w*Dy*Dy
        Dmat(4,ijn) = Dmat(4,ijn) + w*Dy*Dy

        Dmat(6,ijp) = Dmat(6,ijp) + w*Dz*Dz
        Dmat(6,ijn) = Dmat(6,ijn) + w*Dz*Dz  

        Dmat(2,ijp) = Dmat(2,ijp) + w*Dx*Dy
        Dmat(2,ijn) = Dmat(2,ijn) + w*Dx*Dy

        Dmat(3,ijp) = Dmat(3,ijp) + w*Dx*Dz
        Dmat(3,ijn) = Dmat(3,ijn) + w*Dx*Dz
   
        Dmat(5,ijp) = Dmat(5,ijp) + w*Dy*Dz
        Dmat(5,ijn) = Dmat(5,ijn) + w*Dy*Dz     
                                                           
  enddo     

  ! Faces along O-C grid cuts
  do i=1,noc
    ijp = ijl(i)
    ijn = ijr(i)

        w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

        Dx = xc(ijn)-xc(ijp)
        Dy = yc(ijn)-yc(ijp)
        Dz = zc(ijp)-zc(ijn)

        Dmat(1,ijp) = Dmat(1,ijp) + w*Dx*Dx
        Dmat(1,ijn) = Dmat(1,ijn) + w*Dx*Dx 

        Dmat(4,ijp) = Dmat(4,ijp) + w*Dy*Dy
        Dmat(4,ijn) = Dmat(4,ijn) + w*Dy*Dy

        Dmat(6,ijp) = Dmat(6,ijp) + w*Dz*Dz
        Dmat(6,ijn) = Dmat(6,ijn) + w*Dz*Dz  

        Dmat(2,ijp) = Dmat(2,ijp) + w*Dx*Dy
        Dmat(2,ijn) = Dmat(2,ijn) + w*Dx*Dy

        Dmat(3,ijp) = Dmat(3,ijp) + w*Dx*Dz
        Dmat(3,ijn) = Dmat(3,ijn) + w*Dx*Dz
   
        Dmat(5,ijp) = Dmat(5,ijp) + w*Dy*Dz
        Dmat(5,ijn) = Dmat(5,ijn) + w*Dy*Dz 

  end do

  ! Faces on processor boundaries                                             
  do i=1,npro      
      iface = iProcFacesStart + i
      ijp = owner( iface )
      ijn = iProcStart + i

      w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

      Dx = xc(ijn)-xc(ijp)
      Dy = yc(ijn)-yc(ijp)
      Dz = zc(ijn)-zc(ijp)

      Dmat(1,ijp) = Dmat(1,ijp) + w*Dx*Dx

      Dmat(4,ijp) = Dmat(4,ijp) + w*Dy*Dy 

      Dmat(6,ijp) = Dmat(6,ijp) + w*Dz*Dz 

      Dmat(2,ijp) = Dmat(2,ijp) + w*Dx*Dy 

      Dmat(3,ijp) = Dmat(3,ijp) + w*Dx*Dz   
 
      Dmat(5,ijp) = Dmat(5,ijp) + w*Dy*Dz  

  enddo 

  ! Boundary faces:

  do i=numInnerFaces+1,numFaces
      ijp = owner(i)

      w = 1.0_dp/((xf(i)-xc(ijp))**2+(yf(i)-yc(ijp))**2+(zf(i)-zc(ijp))**2)

      Dmat(1,ijp) = Dmat(1,ijp) + w*(xf(i)-xc(ijp))**2
      Dmat(4,ijp) = Dmat(4,ijp) + w*(yf(i)-yc(ijp))**2
      Dmat(6,ijp) = Dmat(6,ijp) + w*(zf(i)-zc(ijp))**2
      Dmat(2,ijp) = Dmat(2,ijp) + w*(xf(i)-xc(ijp))*(yf(i)-yc(ijp)) 
      Dmat(3,ijp) = Dmat(3,ijp) + w*(xf(i)-xc(ijp))*(zf(i)-zc(ijp)) 
      Dmat(5,ijp) = Dmat(5,ijp) + w*(yf(i)-yc(ijp))*(zf(i)-zc(ijp)) 
  end do

  ! Prepare for storage:
  do inp=1,numCells 

    ! Copy from Coefficient matrix 
    D11 = Dmat(1,inp)
    D12 = Dmat(2,inp)
    D13 = Dmat(3,inp)

    D22 = Dmat(4,inp)
    D23 = Dmat(5,inp)
    D33 = Dmat(6,inp)

    ! Symmetric part
    D21 = D12
    D31 = D13
    D32 = D23

    ! Denominator used troughout
    tmp = 1./(d11*d22*d33 - d11*d23*d32 - d12*d21*d33 + d12*d23*d31 + d13*d21*d32 - d13*d22*d31)

    Dmat(1,inp) = (d22*d33 - d23*d32) * tmp
    Dmat(2,inp) = (d21*d33 - d23*d31) * tmp
    Dmat(3,inp) = (d21*d32 - d22*d31) * tmp

    Dmat(4,inp) = (d11*d33 - d13*d31) * tmp
    Dmat(5,inp) = (d12*d33 - d13*d32) * tmp
    Dmat(6,inp) = (d11*d32 - d12*d31) * tmp

    Dmat(7,inp) = (d12*d23 - d13*d22) * tmp
    Dmat(8,inp) = (d11*d23 - d13*d21) * tmp
    Dmat(9,inp) = (d11*d22 - d12*d21) * tmp

  enddo  


!**************************************************************************************************
  elseif(istage.eq.2) then
!**************************************************************************************************

  ! Initialize rhs vector
  b1 = 0.0d0
  b2 = 0.0d0
  b3 = 0.0d0

  ! Inner faces:

  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

    Dx = w * ( xc(ijn)-xc(ijp) ) * ( Fi(ijn)-Fi(ijp) )
    Dy = w * ( yc(ijn)-yc(ijp) ) * ( Fi(ijn)-Fi(ijp) )
    Dz = w * ( zc(ijn)-zc(ijp) ) * ( Fi(ijn)-Fi(ijp) )

    b1(ijp) = b1(ijp) + Dx 
    b1(ijn) = b1(ijn) + Dx

    b2(ijp) = b2(ijp) + Dy
    b2(ijn) = b2(ijn) + Dy 

    b3(ijp) = b3(ijp) + Dz 
    b3(ijn) = b3(ijn) + Dz  
                                                            
  enddo     
  
  ! Faces along O-C grid cuts
  do i=1,noc
    ijp = ijl(i)
    ijn = ijr(i)

    w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

    Dx = w * ( xc(ijn)-xc(ijp) ) * ( Fi(ijn)-Fi(ijp) )
    Dy = w * ( yc(ijn)-yc(ijp) ) * ( Fi(ijn)-Fi(ijp) )
    Dz = w * ( zc(ijn)-zc(ijp) ) * ( Fi(ijn)-Fi(ijp) )

    b1(ijp) = b1(ijp) + Dx 
    b1(ijn) = b1(ijn) + Dx

    b2(ijp) = b2(ijp) + Dy
    b2(ijn) = b2(ijn) + Dy 

    b3(ijp) = b3(ijp) + Dz 
    b3(ijn) = b3(ijn) + Dz
  
  end do


  ! Faces on processor boundaries                                             
  do i=1,npro      
    iface = iProcFacesStart + i
    ijp = owner( iface )
    ijn = iProcStart + i

    w = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)

    Dx = w * ( xc(ijn)-xc(ijp) ) * ( Fi(ijn)-Fi(ijp) )
    Dy = w * ( yc(ijn)-yc(ijp) ) * ( Fi(ijn)-Fi(ijp) )
    Dz = w * ( zc(ijn)-zc(ijp) ) * ( Fi(ijn)-Fi(ijp) )

    b1(ijp) = b1(ijp) + Dx
    b2(ijp) = b2(ijp) + Dy
    b3(ijp) = b3(ijp) + Dz
                                                                                                                        
  enddo

  ! ! Boundary faces:
  ! do i=1,numBoundaryFaces
  !   iface = numInnerFaces + i
  !   ijp = owner(iface)
  !   ijn = numCells+i
  !       w = 1./((xf(i)-xc(ijp))**2+(yf(i)-yc(ijp))**2+(zf(i)-zc(ijp))**2)
  !       b1(ijp) = b1(ijp) + w*(Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
  !       b2(ijp) = b2(ijp) + w*(Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))  
  !       b3(ijp) = b3(ijp) + w*(Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 
  ! end do


  ! Boundary faces:

  do i = 1,ninl
    iface = iInletFacesStart+i
    ijp = owner(iface)
    ijn = iInletStart + i

    w = 1./((xf(i)-xc(ijp))**2+(yf(i)-yc(ijp))**2+(zf(i)-zc(ijp))**2)

    Dx = w*(Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
    Dy = w*(Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))
    Dz = w*(Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 

    b1(ijp) = b1(ijp) + Dx
    b2(ijp) = b2(ijp) + Dy
    b3(ijp) = b3(ijp) + Dz

  enddo

  do i = 1,nout
    iface = iOutletFacesStart+i
    ijp = owner(iface)
    ijn = iOutletStart + i

    w = 1./((xf(i)-xc(ijp))**2+(yf(i)-yc(ijp))**2+(zf(i)-zc(ijp))**2)

    Dx = w*(Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
    Dy = w*(Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))
    Dz = w*(Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 

    b1(ijp) = b1(ijp) + Dx
    b2(ijp) = b2(ijp) + Dy
    b3(ijp) = b3(ijp) + Dz

  enddo

  do i = 1,nsym
    iface = iSymmetryFacesStart+i
    ijp = owner(iface)
    ijn = iSymmetryStart+i

    w = 1./((xf(i)-xc(ijp))**2+(yf(i)-yc(ijp))**2+(zf(i)-zc(ijp))**2)

    Dx = w*(Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
    Dy = w*(Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))
    Dz = w*(Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 

    b1(ijp) = b1(ijp) + Dx
    b2(ijp) = b2(ijp) + Dy
    b3(ijp) = b3(ijp) + Dz

  enddo

  do i = 1,nwal
    iface = iWallFacesStart+i
    ijp = owner(iface)
    ijn = iWallStart+i

    w = 1./((xf(i)-xc(ijp))**2+(yf(i)-yc(ijp))**2+(zf(i)-zc(ijp))**2)

    Dx = w*(Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
    Dy = w*(Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))
    Dz = w*(Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 

    b1(ijp) = b1(ijp) + Dx
    b2(ijp) = b2(ijp) + Dy
    b3(ijp) = b3(ijp) + Dz

  enddo

  do i=1,npru
    iface = iPressOutletFacesStart + i
    ijp = owner(iface)
    ijn = iPressOutletStart + i

    w = 1./((xf(i)-xc(ijp))**2+(yf(i)-yc(ijp))**2+(zf(i)-zc(ijp))**2)

    Dx = w*(Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
    Dy = w*(Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))
    Dz = w*(Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 

    b1(ijp) = b1(ijp) + Dx
    b2(ijp) = b2(ijp) + Dy
    b3(ijp) = b3(ijp) + Dz
    
  enddo


  ! Calculate gradient

  do inp=1,numCells 

    dFidxi(1,inp) = b1(inp)*Dmat(1,inp) - b2(inp)*Dmat(2,inp) +  b3(inp)*Dmat(3,inp) 
    dFidxi(2,inp) = b1(inp)*Dmat(4,inp) - b2(inp)*Dmat(5,inp) -  b3(inp)*Dmat(6,inp) 
    dFidxi(3,inp) = b1(inp)*Dmat(7,inp) - b2(inp)*Dmat(8,inp) +  b3(inp)*Dmat(9,inp) 

  enddo 

   
  endif 
  
end subroutine

