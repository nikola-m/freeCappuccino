subroutine grad_lsq_dm(fi,dFidxi,istage,dmat)
!
!***********************************************************************
!
!      Purpose:
!      Calculates cell-centered gradients using WEIGHTED Least-Squares approach.
!
!      Description:
!      Approach taken from a paper:
!      Dimitry Mavriplis, "Revisiting the Least Squares procedure for Gradient Reconstruction on Unstructured Meshes." NASA/CR-2003-212683.
!      Weights based on inverse cell-centers distance is added to improve conditioning of the system matrix for skewed meshes.
!      Reduces storage requirements compared to QR subroutine.
!      System matrix is symmetric and can be solved efficiently using Cholesky decomposition or by matrix inversion. 
!
!      Arguments:
!
!      FI - field variable which gradient we look for.
!      DFIDXi- cell centered gradient - a three component gradient vector.
!      ISTAGE - integer. If ISTAGE=1 calculates and stores only geometrical
!      parameters - a system matrix for least square problem at every cell. 
!      Usually it is called with ISTAGE=1 at the beggining of simulation.
!      If 2 it doesn't calculate system matrix, just RHS and solves system.
!      Dmat - LSQ matrix with geometry data
!
!      Example call:
!      CALL GRADFI_LSQ_DM(U,DUDX,DUDY,DUDZ,2,D)
!
!***********************************************************************
!
  use types
  use parameters
  use geometry

  implicit none

  integer, intent(in) :: istage
  real(dp),dimension(numTotal), intent(in)   :: fi
  real(dp),dimension(3,numCells), intent(inout) :: dFidxi
  real(dp),dimension(6,numCells), intent(inout) :: dmat

  !
  ! Locals
  !
  integer :: i,ijp,ijn,inp,iface

  real(dp) :: we,ww
  ! real(dp) :: d11,d12,d13,d21,d22,d23,d31,d32,d33
  ! real(dp) :: tmp

  ! For LAPACK DGESV routine
  integer :: INFO
  integer :: IPIV( 3 )
  real(dp) :: A( 3, 3 ), B( 3, 1 )

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

        we = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)
        ww = 1./((xc(ijp)-xc(ijn))**2+(yc(ijp)-yc(ijn))**2+(zc(ijp)-zc(ijn))**2)

        Dmat(1,ijp) = Dmat(1,ijp) + we*(xc(ijn)-xc(ijp))**2
        Dmat(1,ijn) = Dmat(1,ijn) + ww*(xc(ijp)-xc(ijn))**2 

        Dmat(4,ijp) = Dmat(4,ijp) + we*(yc(ijn)-yc(ijp))**2
        Dmat(4,ijn) = Dmat(4,ijn) + ww*(yc(ijp)-yc(ijn))**2

        Dmat(6,ijp) = Dmat(6,ijp) + we*(zc(ijn)-zc(ijp))**2
        Dmat(6,ijn) = Dmat(6,ijn) + ww*(zc(ijp)-zc(ijn))**2  

        Dmat(2,ijp) = Dmat(2,ijp) + we*(xc(ijn)-xc(ijp))*(yc(ijn)-yc(ijp)) 
        Dmat(2,ijn) = Dmat(2,ijn) + ww*(xc(ijp)-xc(ijn))*(yc(ijp)-yc(ijn))

        Dmat(3,ijp) = Dmat(3,ijp) + we*(xc(ijn)-xc(ijp))*(zc(ijn)-zc(ijp)) 
        Dmat(3,ijn) = Dmat(3,ijn) + ww*(xc(ijp)-xc(ijn))*(zc(ijp)-zc(ijn))
   
        Dmat(5,ijp) = Dmat(5,ijp) + we*(yc(ijn)-yc(ijp))*(zc(ijn)-zc(ijp)) 
        Dmat(5,ijn) = Dmat(5,ijn) + ww*(yc(ijp)-yc(ijn))*(zc(ijp)-zc(ijn))    
                                                           
  enddo     

  ! Faces along O-C grid cuts
  do i=1,noc
    ijp = ijl(i)
    ijn = ijr(i)

        we = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)
        ww = 1./((xc(ijp)-xc(ijn))**2+(yc(ijp)-yc(ijn))**2+(zc(ijp)-zc(ijn))**2)

        Dmat(1,ijp) = Dmat(1,ijp) + we*(xc(ijn)-xc(ijp))**2
        Dmat(1,ijn) = Dmat(1,ijn) + ww*(xc(ijp)-xc(ijn))**2 

        Dmat(4,ijp) = Dmat(4,ijp) + we*(yc(ijn)-yc(ijp))**2
        Dmat(4,ijn) = Dmat(4,ijn) + ww*(yc(ijp)-yc(ijn))**2

        Dmat(6,ijp) = Dmat(6,ijp) + we*(zc(ijn)-zc(ijp))**2
        Dmat(6,ijn) = Dmat(6,ijn) + ww*(zc(ijp)-zc(ijn))**2  

        Dmat(2,ijp) = Dmat(2,ijp) + we*(xc(ijn)-xc(ijp))*(yc(ijn)-yc(ijp)) 
        Dmat(2,ijn) = Dmat(2,ijn) + ww*(xc(ijp)-xc(ijn))*(yc(ijp)-yc(ijn))

        Dmat(3,ijp) = Dmat(3,ijp) + we*(xc(ijn)-xc(ijp))*(zc(ijn)-zc(ijp)) 
        Dmat(3,ijn) = Dmat(3,ijn) + ww*(xc(ijp)-xc(ijn))*(zc(ijp)-zc(ijn))
   
        Dmat(5,ijp) = Dmat(5,ijp) + we*(yc(ijn)-yc(ijp))*(zc(ijn)-zc(ijp)) 
        Dmat(5,ijn) = Dmat(5,ijn) + ww*(yc(ijp)-yc(ijn))*(zc(ijp)-zc(ijn)) 

  end do

  ! Boundary faces:

  do i=numInnerFaces+1,numFaces
      ijp = owner(i)

      we = 1.0_dp/((xf(i)-xc(ijp))**2+(yf(i)-yc(ijp))**2+(zf(i)-zc(ijp))**2)

      Dmat(1,ijp) = Dmat(1,ijp) + we*(xf(i)-xc(ijp))**2
      Dmat(4,ijp) = Dmat(4,ijp) + we*(yf(i)-yc(ijp))**2
      Dmat(6,ijp) = Dmat(6,ijp) + we*(zf(i)-zc(ijp))**2
      Dmat(2,ijp) = Dmat(2,ijp) + we*(xf(i)-xc(ijp))*(yf(i)-yc(ijp)) 
      Dmat(3,ijp) = Dmat(3,ijp) + we*(xf(i)-xc(ijp))*(zf(i)-zc(ijp)) 
      Dmat(5,ijp) = Dmat(5,ijp) + we*(yf(i)-yc(ijp))*(zf(i)-zc(ijp)) 
  end do


!**************************************************************************************************
  elseif(istage.eq.2) then
!**************************************************************************************************

  ! Initialize rhs vector
  b1(:) = 0.0d0
  b2(:) = 0.0d0
  b3(:) = 0.0d0

  ! Inner faces:

  do i=1,numInnerFaces                                                       
        ijp = owner(i)
        ijn = neighbour(i)

        we = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)
        ww = 1./((xc(ijp)-xc(ijn))**2+(yc(ijp)-yc(ijn))**2+(zc(ijp)-zc(ijn))**2)

        b1(ijp) = b1(ijp) + we*(Fi(ijn)-Fi(ijp))*(xc(ijn)-xc(ijp)) 
        b1(ijn) = b1(ijn) + ww*(Fi(ijp)-Fi(ijn))*(xc(ijp)-xc(ijn))

        b2(ijp) = b2(ijp) + we*(Fi(ijn)-Fi(ijp))*(yc(ijn)-yc(ijp)) 
        b2(ijn) = b2(ijn) + ww*(Fi(ijp)-Fi(ijn))*(yc(ijp)-yc(ijn))  

        b3(ijp) = b3(ijp) + we*(Fi(ijn)-Fi(ijp))*(zc(ijn)-zc(ijp)) 
        b3(ijn) = b3(ijn) + ww*(Fi(ijp)-Fi(ijn))*(zc(ijp)-zc(ijn))
                                                            
  enddo     
  
  ! Faces along O-C grid cuts
  do i=1,noc
    ijp = ijl(i)
    ijn = ijr(i)

        we = 1./((xc(ijn)-xc(ijp))**2+(yc(ijn)-yc(ijp))**2+(zc(ijn)-zc(ijp))**2)
        ww = 1./((xc(ijp)-xc(ijn))**2+(yc(ijp)-yc(ijn))**2+(zc(ijp)-zc(ijn))**2)
        
        b1(ijp) = b1(ijp) + we*(Fi(ijn)-Fi(ijp))*(xc(ijn)-xc(ijp)) 
        b1(ijn) = b1(ijn) + ww*(Fi(ijp)-Fi(ijn))*(xc(ijp)-xc(ijn))

        b2(ijp) = b2(ijp) + we*(Fi(ijn)-Fi(ijp))*(yc(ijn)-yc(ijp)) 
        b2(ijn) = b2(ijn) + ww*(Fi(ijp)-Fi(ijn))*(yc(ijp)-yc(ijn))  

        b3(ijp) = b3(ijp) + we*(Fi(ijn)-Fi(ijp))*(zc(ijn)-zc(ijp)) 
        b3(ijn) = b3(ijn) + ww*(Fi(ijp)-Fi(ijn))*(zc(ijp)-zc(ijn))
  
  end do

  ! Boundary faces:
  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)
    ijn = numCells+i
        we = 1./((xf(i)-xc(ijp))**2+(yf(i)-yc(ijp))**2+(zf(i)-zc(ijp))**2)
        b1(ijp) = b1(ijp) + we*(Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
        b2(ijp) = b2(ijp) + we*(Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))  
        b3(ijp) = b3(ijp) + we*(Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 
  end do


  ! Calculate gradient

  ! Cell loop
  do inp=1,numCells

        ! ! Copy from Coefficient matrix 
        ! d11 = Dmat(1,inp)
        ! d12 = Dmat(2,inp)
        ! d13 = Dmat(3,inp)

        ! d22 = Dmat(4,inp)
        ! d23 = Dmat(5,inp)
        ! d33 = Dmat(6,inp)

        ! ! Symmetric part
        ! d21 = d12
        ! d31 = d13
        ! d32 = d23

        ! ! Solve system 

        ! tmp = 1./(d11*d22*d33 - d11*d23*d32 - d12*d21*d33 + d12*d23*d31 + d13*d21*d32 - d13*d22*d31)

        ! dFidxi(1,inp) = ( (b1(inp)*(d22*d33 - d23*d32)) - (b2(inp)*(d21*d33 - d23*d31)) + (b3(inp)*(d21*d32 - d22*d31)) ) * tmp
        ! dFidxi(2,inp) = ( (b2(inp)*(d11*d33 - d13*d31)) - (b1(inp)*(d12*d33 - d13*d32)) - (b3(inp)*(d11*d32 - d12*d31)) ) * tmp
        ! dFidxi(3,inp) = ( (b1(inp)*(d12*d23 - d13*d22)) - (b2(inp)*(d11*d23 - d13*d21)) + (b3(inp)*(d11*d22 - d12*d21)) ) * tmp

!
!     Solve the system A*X = B.
!
      A(1,1) = Dmat(1,inp)
      A(1,2) = Dmat(2,inp)
      A(1,3) = Dmat(3,inp)
      A(2,1) = A(1,2)
      A(2,2) = Dmat(4,inp)
      A(2,3) = Dmat(5,inp)
      A(3,1) = A(1,3)
      A(3,2) = A(2,3)
      A(3,3) = Dmat(6,inp)

      B(1,1) = b1(inp)
      B(2,1) = b2(inp)
      B(3,1) = b3(inp)

      CALL DGESV( 3, 1, A, 3, IPIV, B, 3, INFO )

      dFidxi(1,inp) = B(1,1)
      dFidxi(2,inp) = B(2,1)
      dFidxi(3,inp) = B(3,1)

  enddo
   
  endif 
  
end subroutine

