! Global
!     continuityErrs
!
! Description
!     Calculates and prints the continuity errors.
!---------------------------------------------------------------------------


! Calculate Volume Scalar Field equal to the divergence of Phi (the mass flow trough facce)

  ! Initialize array with zero value.
  su(:) = 0.0d0

  ! Inner faces                                               
  do i=1,numInnerFaces                                                    
    ijp = owner(i)
    ijn = neighbour(i)
    su(ijp) = su(ijp)-flmass(ijp)
    su(ijn) = su(ijn)+flmass(ijp)
  enddo     
  
  ! Faces along O-C grid cuts
  do i=1,noc
    ijp = ijl(i)
    ijn = ijr(i)
    su(ijp)=su(ijp)-fmoc(i)
    su(ijn)=su(ijn)+fmoc(i)
  end do

  ! Inlet boundaries (mass fluxes prescribed in routine 'bcin')
  do i=1,ninl
    ijp = owner(iInletFacesStart+i)
    su(ijp)=su(ijp)+fmi(i)
  end do

  ! Correct mass flux to satisfy global mass conservation & add to source
  do i=1,nout
    ijp = owner(iOutletFacesStart+i)
    su(ijp)=su(ijp)-fmo(i)
  end do


  sumLocalContErr = sum( abs( su ) ) 

   globalContErr = sum( su )

   cumulativeContErr = cumulativeContErr + globalContErr


  write(6,'(3(a,es10.3))') "  time step continuity errors : sum local = ", sumLocalContErr, &
 &                          ", global = ", globalContErr, &
 &                          ", cumulative = ", cumulativeContErr

