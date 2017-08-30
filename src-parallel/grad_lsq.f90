subroutine grad_lsq(fi,dFidxi,istage,dmat)
!
!***********************************************************************
!
!      Purpose:
!      Calculates cell-centered gradients using UNWEIGHTED Least-Squares approach.
!
!      Description:
!      Approach taken from PhD thesis of Bojan Niceno, TU Delft, 2000.,
!      also in Muzaferija and Gossman JCP paper from 1995.
!
!      Arguments:
!
!      FI - field variable which gradient we look for.
!      DFiDXi - cell centered gradient - a three component gradient vector.
!      ISTAGE - integer. If ISTAGE=1 calculates and stores only geometrical
!      parameters - a system matrix for least square problem at every cell. 
!      Usually it is called with ISTAGE=1 at the beggining of simulation.
!      If 2 it doesn't calculate system matrix, just RHS and solves system.
!      Dmat - LSQ matrix with geometry data
!
!      Example call:
!      CALL GRADFI_LSQ(U,DUDXI,2,D)
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
  real(dp),dimension(6,numCells), intent(inout) :: Dmat

  !
  ! Locals
  !
  integer :: i,ijp,ijn,inp,iface

  ! real(dp), dimension(numCells) :: b1,b2,b3
  real(dp) :: Dx,Dy,Dz,Jac,Dinv(6),DPhi1,DPhi2 ! Dxc2,Dyc2,Dzc2

  ! ! For LAPACK DGESV routine
  ! integer :: INFO
  ! ! integer :: IPIV( 3 )
  ! real(dp) :: A( 3, 3 ), B( 3, 1 )
!
!***********************************************************************
!

  ! Initialize dmat matrix:
  Dmat = 0.0_dp

  if(istage.eq.1) then
  ! Coefficient matrix - should be calculated only once 


  ! Inner faces:                                             
  do i=1,numInnerFaces                                                       
      ijp = owner(i)
      ijn = neighbour(i)

      Dx = xc(ijn)-xc(ijp)
      Dy = yc(ijn)-yc(ijp)
      Dz = zc(ijn)-zc(ijp)

      ! Dxc2 = xc(ijn)-xc(ijp)
      ! Dyc2 = yc(ijn)-yc(ijp)
      ! Dzc2 = zc(ijn)-zc(ijp)

      Dmat(1,ijp) = Dmat(1,ijp) + Dx*Dx  ! 1,1
      Dmat(1,ijn) = Dmat(1,ijn) + Dx*Dx

      Dmat(2,ijp) = Dmat(2,ijp) + Dy*Dy  ! 2,2 
      Dmat(2,ijn) = Dmat(2,ijn) + Dy*Dy 

      Dmat(3,ijp) = Dmat(3,ijp) + Dz*Dz  ! 3,3 
      Dmat(3,ijn) = Dmat(3,ijn) + Dz*Dz  

      Dmat(4,ijp) = Dmat(4,ijp) + Dx*Dy  ! 1,2  &  2,1 
      Dmat(4,ijn) = Dmat(4,ijn) + Dx*Dy  

      Dmat(5,ijp) = Dmat(5,ijp) + Dx*Dz  ! 1,3  &  3,1  
      Dmat(5,ijn) = Dmat(5,ijn) + Dx*Dz 
 
      Dmat(6,ijp) = Dmat(6,ijp) + Dy*Dz  ! 2,3  &  3,2 
      Dmat(6,ijn) = Dmat(6,ijn) + Dy*Dz                                                                  
enddo     

  ! ! Inner faces:                                             
  ! do i=1,numInnerFaces                                                       
  !       ijp = owner(i)
  !       ijn = neighbour(i)

  !       Dmat(1,ijp) = Dmat(1,ijp) + (xc(ijn)-xc(ijp))**2
  !       Dmat(1,ijn) = Dmat(1,ijn) + (xc(ijp)-xc(ijn))**2 

  !       Dmat(4,ijp) = Dmat(4,ijp) + (yc(ijn)-yc(ijp))**2
  !       Dmat(4,ijn) = Dmat(4,ijn) + (yc(ijp)-yc(ijn))**2

  !       Dmat(6,ijp) = Dmat(6,ijp) + (zc(ijn)-zc(ijp))**2
  !       Dmat(6,ijn) = Dmat(6,ijn) + (zc(ijp)-zc(ijn))**2  

  !       Dmat(2,ijp) = Dmat(2,ijp) + (xc(ijn)-xc(ijp))*(yc(ijn)-yc(ijp)) 
  !       Dmat(2,ijn) = Dmat(2,ijn) + (xc(ijp)-xc(ijn))*(yc(ijp)-yc(ijn))

  !       Dmat(3,ijp) = Dmat(3,ijp) + (xc(ijn)-xc(ijp))*(zc(ijn)-zc(ijp)) 
  !       Dmat(3,ijn) = Dmat(3,ijn) + (xc(ijp)-xc(ijn))*(zc(ijp)-zc(ijn))
   
  !       Dmat(5,ijp) = Dmat(5,ijp) + (yc(ijn)-yc(ijp))*(zc(ijn)-zc(ijp)) 
  !       Dmat(5,ijn) = Dmat(5,ijn) + (yc(ijp)-yc(ijn))*(zc(ijp)-zc(ijn))                                                                 
  ! enddo     
 
  
  ! ! Faces along O-C grid cuts
  ! do i=1,noc
  !   ijp = ijl(i)
  !   ijn = ijr(i)

  !       Dmat(1,ijp) = Dmat(1,ijp) + (xc(ijn)-xc(ijp))**2
  !       Dmat(1,ijn) = Dmat(1,ijn) + (xc(ijp)-xc(ijn))**2 

  !       Dmat(4,ijp) = Dmat(4,ijp) + (yc(ijn)-yc(ijp))**2
  !       Dmat(4,ijn) = Dmat(4,ijn) + (yc(ijp)-yc(ijn))**2

  !       Dmat(6,ijp) = Dmat(6,ijp) + (zc(ijn)-zc(ijp))**2
  !       Dmat(6,ijn) = Dmat(6,ijn) + (zc(ijp)-zc(ijn))**2  

  !       Dmat(2,ijp) = Dmat(2,ijp) + (xc(ijn)-xc(ijp))*(yc(ijn)-yc(ijp)) 
  !       Dmat(2,ijn) = Dmat(2,ijn) + (xc(ijp)-xc(ijn))*(yc(ijp)-yc(ijn))

  !       Dmat(3,ijp) = Dmat(3,ijp) + (xc(ijn)-xc(ijp))*(zc(ijn)-zc(ijp)) 
  !       Dmat(3,ijn) = Dmat(3,ijn) + (xc(ijp)-xc(ijn))*(zc(ijp)-zc(ijn))
   
  !       Dmat(5,ijp) = Dmat(5,ijp) + (yc(ijn)-yc(ijp))*(zc(ijn)-zc(ijp)) 
  !       Dmat(5,ijn) = Dmat(5,ijn) + (yc(ijp)-yc(ijn))*(zc(ijp)-zc(ijn))  

  ! end do

  
  ! Boundary faces:

  do i=1,numBoundaryFaces
    iface = numInnerFaces + i
    ijp = owner(iface)

        Dx = xf(iface)-xc(ijp)
        Dy = yf(iface)-yc(ijp)
        Dz = zf(iface)-zc(ijp)

        Dmat(1,ijp) = Dmat(1,ijp) + Dx*Dx
        Dmat(2,ijp) = Dmat(2,ijp) + Dy*Dy
        Dmat(3,ijp) = Dmat(3,ijp) + Dz*Dz
        Dmat(4,ijp) = Dmat(4,ijp) + Dx*Dy 
        Dmat(5,ijp) = Dmat(5,ijp) + Dx*Dz 
        Dmat(6,ijp) = Dmat(6,ijp) + Dy*Dz 
  end do

  ! ! Inlet: 
  ! do i=1,ninl
  !   iface = iInletFacesStart + i
  !   ijp = owner(iface)
  !       Dmat(1,ijp) = Dmat(1,ijp) + (xf(iface)-xc(ijp))**2
  !       Dmat(4,ijp) = Dmat(4,ijp) + (yf(iface)-yc(ijp))**2
  !       Dmat(6,ijp) = Dmat(6,ijp) + (zf(iface)-zc(ijp))**2
  !       Dmat(2,ijp) = Dmat(2,ijp) + (xf(iface)-xc(ijp))*(yf(iface)-yc(ijp)) 
  !       Dmat(3,ijp) = Dmat(3,ijp) + (xf(iface)-xc(ijp))*(zf(iface)-zc(ijp)) 
  !       Dmat(5,ijp) = Dmat(5,ijp) + (yf(iface)-yc(ijp))*(zf(iface)-zc(ijp)) 
  ! end do

  ! ! Outlet
  ! do i=1,nout
  ! iface = iOutletFacesStart + i
  ! ijp = owner(iface)
  !       Dmat(1,ijp) = Dmat(1,ijp) + (xf(iface)-xc(ijp))**2
  !       Dmat(4,ijp) = Dmat(4,ijp) + (yf(iface)-yc(ijp))**2
  !       Dmat(6,ijp) = Dmat(6,ijp) + (zf(iface)-zc(ijp))**2
  !       Dmat(2,ijp) = Dmat(2,ijp) + (xf(iface)-xc(ijp))*(yf(iface)-yc(ijp)) 
  !       Dmat(3,ijp) = Dmat(3,ijp) + (xf(iface)-xc(ijp))*(zf(iface)-zc(ijp)) 
  !       Dmat(5,ijp) = Dmat(5,ijp) + (yf(iface)-yc(ijp))*(zf(iface)-zc(ijp)) 
  ! end do

  ! ! Symmetry
  ! do i=1,nsym
  ! iface = iSymmetryFacesStart + i
  ! ijp = owner(iface)
  !       Dmat(1,ijp) = Dmat(1,ijp) + (xf(iface)-xc(ijp))**2
  !       Dmat(4,ijp) = Dmat(4,ijp) + (yf(iface)-yc(ijp))**2
  !       Dmat(6,ijp) = Dmat(6,ijp) + (zf(iface)-zc(ijp))**2
  !       Dmat(2,ijp) = Dmat(2,ijp) + (xf(iface)-xc(ijp))*(yf(iface)-yc(ijp)) 
  !       Dmat(3,ijp) = Dmat(3,ijp) + (xf(iface)-xc(ijp))*(zf(iface)-zc(ijp)) 
  !       Dmat(5,ijp) = Dmat(5,ijp) + (yf(iface)-yc(ijp))*(zf(iface)-zc(ijp)) 
  ! end do

  ! ! Wall
  ! do i=1,nwal
  ! iface = iWallFacesStart + i
  ! ijp = owner(iface)
  !       Dmat(1,ijp) = Dmat(1,ijp) + (xf(iface)-xc(ijp))**2
  !       Dmat(4,ijp) = Dmat(4,ijp) + (yf(iface)-yc(ijp))**2
  !       Dmat(6,ijp) = Dmat(6,ijp) + (zf(iface)-zc(ijp))**2
  !       Dmat(2,ijp) = Dmat(2,ijp) + (xf(iface)-xc(ijp))*(yf(iface)-yc(ijp)) 
  !       Dmat(3,ijp) = Dmat(3,ijp) + (xf(iface)-xc(ijp))*(zf(iface)-zc(ijp)) 
  !       Dmat(5,ijp) = Dmat(5,ijp) + (yf(iface)-yc(ijp))*(zf(iface)-zc(ijp)) 
  ! end do

  ! ! Pressure Outlet
  ! do i=1,npru
  ! iface = iPressOutletFacesStart + i
  ! ijp = owner(iface)
  !       Dmat(1,ijp) = Dmat(1,ijp) + (xf(iface)-xc(ijp))**2
  !       Dmat(4,ijp) = Dmat(4,ijp) + (yf(iface)-yc(ijp))**2
  !       Dmat(6,ijp) = Dmat(6,ijp) + (zf(iface)-zc(ijp))**2
  !       Dmat(2,ijp) = Dmat(2,ijp) + (xf(iface)-xc(ijp))*(yf(iface)-yc(ijp)) 
  !       Dmat(3,ijp) = Dmat(3,ijp) + (xf(iface)-xc(ijp))*(zf(iface)-zc(ijp)) 
  !       Dmat(5,ijp) = Dmat(5,ijp) + (yf(iface)-yc(ijp))*(zf(iface)-zc(ijp)) 
  ! end do

  ! Find and store the inverse of LSQ matrix Dmat
  do inp=1,numCells 

    Jac  =   Dmat(1,inp) * Dmat(2,inp) * Dmat(3,inp)                         &
     -       Dmat(1,inp) * Dmat(6,inp) * Dmat(6,inp)                         &
     -       Dmat(4,inp) * Dmat(4,inp) * Dmat(3,inp)                         &
     +   2 * Dmat(4,inp) * Dmat(5,inp) * Dmat(6,inp)                         &
     -       Dmat(5,inp) * Dmat(5,inp) * Dmat(2,inp)

    Dinv(1) = +( Dmat(2,inp)*Dmat(3,inp) - Dmat(6,inp)*Dmat(6,inp) ) / (Jac+small)
    Dinv(2) = +( Dmat(1,inp)*Dmat(3,inp) - Dmat(5,inp)*Dmat(5,inp) ) / (Jac+small)
    Dinv(3) = +( Dmat(1,inp)*Dmat(2,inp) - Dmat(4,inp)*Dmat(4,inp) ) / (Jac+small)
    Dinv(4) = -( Dmat(4,inp)*Dmat(3,inp) - Dmat(5,inp)*Dmat(6,inp) ) / (Jac+small)
    Dinv(5) = +( Dmat(4,inp)*Dmat(6,inp) - Dmat(5,inp)*Dmat(2,inp) ) / (Jac+small)
    Dinv(6) = -( Dmat(1,inp)*Dmat(6,inp) - Dmat(4,inp)*Dmat(5,inp) ) / (Jac+small)

    Dmat(1,inp) = Dinv(1) 
    Dmat(2,inp) = Dinv(2)
    Dmat(3,inp) = Dinv(3)
    Dmat(4,inp) = Dinv(4)
    Dmat(5,inp) = Dinv(5)
    Dmat(6,inp) = Dinv(6)

  end do 

  elseif(istage.eq.2) then

  ! ! Initialize rhs vector
  ! b1 = 0.0_dp
  ! b2 = 0.0_dp
  ! b3 = 0.0_dp




  ! Inner faces:

  do i=1,numInnerFaces                                                       
    ijp = owner(i)
    ijn = neighbour(i)

    Dx = xc(ijn)-xc(ijp)
    Dy = yc(ijn)-yc(ijp)
    Dz = zc(ijn)-zc(ijp)

    DPhi1 = Fi(ijn)-Fi(ijp)
    DPhi2 = Fi(ijn)-Fi(ijp)

    dFidxi(1,ijp) = dFidxi(1,ijp) + DPhi1*(Dmat(1,ijp)*Dx+Dmat(4,ijp)*Dy+Dmat(5,ijp)*Dz)
    dFidxi(1,ijn) = dFidxi(1,ijn) + DPhi2*(Dmat(1,ijn)*Dx+Dmat(4,ijn)*Dy+Dmat(5,ijn)*Dz)

    dFidxi(2,ijp) = dFidxi(2,ijp) + DPhi1*(Dmat(4,ijp)*Dx+Dmat(2,ijp)*Dy+Dmat(6,ijp)*Dz)
    dFidxi(2,ijn) = dFidxi(2,ijn) + DPhi2*(Dmat(4,ijn)*Dx+Dmat(2,ijn)*Dy+Dmat(6,ijn)*Dz)

    dFidxi(3,ijp) = dFidxi(3,ijp) + DPhi1*(Dmat(5,ijp)*Dx+Dmat(6,ijp)*Dy+Dmat(3,ijp)*Dz)
    dFidxi(3,ijn) = dFidxi(3,ijn) + DPhi2*(Dmat(5,ijn)*Dx+Dmat(6,ijn)*Dy+Dmat(3,ijn)*Dz)                                                                                                                       
  enddo     

  ! do i=1,numInnerFaces                                                       
  !       ijp = owner(i)
  !       ijn = neighbour(i)

  !       b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*Dx 
  !       b1(ijn) = b1(ijn) + (Fi(ijp)-Fi(ijn))*Dx

  !       b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*Dy 
  !       b2(ijn) = b2(ijn) + (Fi(ijp)-Fi(ijn))*Dy 

  !       b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*Dz
  !       b3(ijn) = b3(ijn) + (Fi(ijp)-Fi(ijn))*Dz                                                                                                                       
  ! enddo     


! Kod Nicena, ako je u pitanju granica i ako je simetrija imamo nulu, 
! prirodno jer su Fi isti sa obe strane

  
  ! ! Faces along O-C grid cuts
  ! do i=1,noc
  !   ijp = ijl(i)
  !   ijn = ijr(i)

  !       b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xc(ijn)-xc(ijp)) 
  !       b1(ijn) = b1(ijn) + (Fi(ijp)-Fi(ijn))*(xc(ijp)-xc(ijn))

  !       b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yc(ijn)-yc(ijp)) 
  !       b2(ijn) = b2(ijn) + (Fi(ijp)-Fi(ijn))*(yc(ijp)-yc(ijn))  

  !       b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zc(ijn)-zc(ijp)) 
  !       b3(ijn) = b3(ijn) + (Fi(ijp)-Fi(ijn))*(zc(ijp)-zc(ijn))
  
  ! end do

  ! Boundary faces:

  ! do i = 1,ninl
  !   iface = iInletFacesStart+i
  !   ijp = owner(iface)
  !   ijn = iInletStart + i
  !       b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
  !       b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))  
  !       b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 
  ! enddo

  ! do i = 1,nout
  !   iface = iOutletFacesStart+i
  !   ijp = owner(iface)
  !   ijn = iOutletStart + i
  !       b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
  !       b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))  
  !       b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 
  ! enddo

  ! do i = 1,nsym
  !   iface = iSymmetryFacesStart+i
  !   ijp = owner(iface)
  !   ijn = iSymmetryStart+i
  !       b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
  !       b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))  
  !       b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 
  ! enddo

  ! do i = 1,nwal
  !   iface = iWallFacesStart+i
  !   ijp = owner(iface)
  !   ijn = iWallStart+i
  !       b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
  !       b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))  
  !       b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 
  ! enddo

  ! do i=1,npru
  !   iface = iPressOutletFacesStart + i
  !   ijp = owner(iface)
  !   ijn = iPressOutletStart + i
  !       b1(ijp) = b1(ijp) + (Fi(ijn)-Fi(ijp))*(xf(iface)-xc(ijp)) 
  !       b2(ijp) = b2(ijp) + (Fi(ijn)-Fi(ijp))*(yf(iface)-yc(ijp))  
  !       b3(ijp) = b3(ijp) + (Fi(ijn)-Fi(ijp))*(zf(iface)-zc(ijp)) 
  ! enddo


  ! !
  ! ! Solve the system A*X = B.
  ! !

  ! do inp=1,numCells 

  !   A(1,1) = Dmat(1,inp)
  !   A(1,2) = Dmat(2,inp)
  !   A(1,3) = Dmat(3,inp)
  !   A(2,2) = Dmat(4,inp)
  !   A(2,3) = Dmat(5,inp)
  !   A(3,3) = Dmat(6,inp)

  !   ! ! Only if we use DGESV linear solver:
  !   ! A(2,1) = A(1,2)
  !   ! A(3,1) = A(1,3)
  !   ! A(3,2) = A(2,3)


  !   B(1,1) = b1(inp)
  !   B(2,1) = b2(inp)
  !   B(3,1) = b3(inp)

  !   ! Solve 3x3 linear system:
  !   ! CALL DGESV( 3, 1, A, 3, IPIV, B, 3, INFO )

  !   ! ... or exploit symmetry:
  !   CALL DPOSV( 'Upper', 3, 1, A, 3, B, 3, INFO )

  !   dFidxi(1,inp) = B(1,1)
  !   dFidxi(2,inp) = B(2,1)
  !   dFidxi(3,inp) = B(3,1)

  ! enddo
   

endif 
  
end subroutine