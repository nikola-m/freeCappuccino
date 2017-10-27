! Global
!     continuityErrs
!
! Description
!     Calculates and prints the continuity errors.
!---------------------------------------------------------------------------


! Calculate Volume Scalar Field equal to the divergence of Phi (the mass flow trough facce)

  ! Initialize array with zero value.
  res = 0.0_dp

  ! Inner faces                                               
  do i=1,numInnerFaces                                                    
    ijp = owner(i)
    ijn = neighbour(i)
    res(ijp) = res(ijp)-flmass(ijp)
    res(ijn) = res(ijn)+flmass(ijp)
  enddo  
  
  ! Faces along O-C grid cuts
  do i=1,noc
    ijp = ijl(i)
    ijn = ijr(i)
    res(ijp)=res(ijp)-fmoc(i)
    res(ijn)=res(ijn)+fmoc(i)
  end do


  ! Inlet boundaries (mass fluxes prescribed in routine 'bcin')
  do i=1,ninl
    ijp = owner(iInletFacesStart+i)
    res(ijp)=res(ijp)-fmi(i)
  end do


  ! Correct mass flux to satisfy global mass conservation
  do i=1,nout
    ijp = owner(iOutletFacesStart+i)
    res(ijp)=res(ijp)-fmo(i)
  end do

  ! The way it is done in OpenFOAM:
  ! sumLocalContErr = timestep*volumeWeightedAverage( abs(res) )
  ! globalContErr = timestep*volumeWeightedAverage( res )

  sumLocalContErr = sum( abs( res ) ) 
  
  globalContErr = sum( res )

  cumulativeContErr = cumulativeContErr + globalContErr


  write(6,'(3(a,es10.3))') "  time step continuity errors : sum local = ", sumLocalContErr, &
 &                          ", global = ", globalContErr, &
 &                          ", cumulative = ", cumulativeContErr

