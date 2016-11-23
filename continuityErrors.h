! Global
!     continuityErrs
!
! Description
!     Calculates and prints the continuity errors.
!---------------------------------------------------------------------------


!//     volScalarField contErr(fvc::div(phi));


! Calculate Volume Scalar Field equal to the divergence of Phi (the mass flow trough facce)

  ! Initialize array with zero value.
  su(:) = 0.0d0

  ! Inner faces                                               
  do i=1,numInnerFaces                                                    
    ijp = owner(i)
    ijn = neighbour(i)
    are = sqrt( arx(i)**2 + ary(i)**2 + arz(i)**2 )
    su(ijp) = su(ijp)-flmass(ijp)*are
    su(ijn) = su(ijn)+flmass(ijp)*are                                                                                                                           
  enddo     
  
  ! Faces along O-C grid cuts
  do i=1,noc
    ijp = ijl(i)
    ijn = ijr(i)
    are = sqrt( xnoc(i)**2 + ynoc(i)**2 +znoc(i)**2 )
    su(ijp)=su(ijp)-fmoc(i)*are
    su(ijn)=su(ijn)+fmoc(i)*are
  end do

  ! Inlet boundaries (mass fluxes prescribed in routine 'bcin')
  do i=1,ninl
    ijp = owner(iInletFacesStart+i)
    are = sqrt( xni(i)**2 + yni(i)**2 +zni(i)**2 )
    su(ijp)=su(ijp)+fmi(i)*are
  end do

  ! Correct mass flux to satisfy global mass conservation & add to source
  do i=1,nout
    ijp = owner(iOutletFacesStart+i)
    are = sqrt( xno(i)**2 + yno(i)**2 +zno(i)**2 )
    su(ijp)=su(ijp)-fmo(i)*are
  end do


!//     scalar sumLocalContErr = runTime.deltaTValue()*
!//         mag(contErr)().weightedAverage(mesh.V()).value();
      sumLocalContErr = timestep*volumeWeightedAverage(abs(su))


!//     scalar globalContErr = runTime.deltaTValue()*
!//         contErr.weightedAverage(mesh.V()).value();
       globalContErr = timestep*volumeWeightedAverage(su)

!//     cumulativeContErr += globalContErr;
       cumulativeContErr = cumulativeContErr + globalContErr


      write(66,'(3(a,es10.3))') "time step continuity errors : sum local = ", sumLocalContErr, &
     &                          ", global = ", globalContErr, &
     &                          ", cumulative = ", cumulativeContErr

!// ************************************************************************* //

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

  ! Outlet boundaries
  do i=1,nout
    ijp = owner(iOutletFacesStart+i)
    su(ijp)=su(ijp)-fmo(i)
  end do

  write(66,'(20x,a,1pe10.3,/,20x,a,1pe10.3)') ' sum  =',sum( su(:) ), &
                                              '|sum| =',sum( abs(su(:)) )
