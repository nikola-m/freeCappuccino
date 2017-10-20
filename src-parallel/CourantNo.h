!
! Calculate and output the mean and maximum Courant Numbers.
!

 if (ltransient) then
 CoNum = 0.0_dp
 meanCoNum = 0.0_dp

 su = 0.0_dp

  !
  ! Suface sum of magnitude (i.e. absolute value) of mass flux phi, over inner faces only (which includes o- c- faces)
  !

  ! Inner faces                                               
  do i=1,numInnerFaces                                                                                                               
    ijp = owner(i)
    ijn = neighbour(i)
    su(ijp) = su(ijp)+abs(flmass(i))
    su(ijn) = su(ijn)+abs(flmass(i))                                                                                                               
  enddo     

  
  ! Faces along O-C grid cuts
  do i=1,noc
    ijp = ijl(i)
    ijn = ijr(i)
    su(ijp) = su(ijp)+abs(fmoc(i))
    su(ijn) = su(ijn)+abs(fmoc(i))
  end do

  ! Faces on processor boundaries
  do i=1,npro
    ijp = owner( iProcFacesStart + i )
    su(ijp) = su(ijp)+abs(fmpro(i))
  end do

  !// Mozda ne idu jer ide samo internal field...
  !// ! Inlet boundaries (mass fluxes prescribed in routine 'bcin')
  !// do i=1,ninl
  !//   ijp = owner(iInletFacesStart+i)
  !//   su(ijp)=su(ijp)+abs(fmi(i))
  !// end do
  !// ! Outlet boundaries
  !// do i=1,nout
  !//   ijp = owner(iOutletFacesStart+i)
  !//   su(ijp)=su(ijp)+abs(fmo(i))
  !// end do

  ! Accumulate by looping trough cells
  suma = 0.0_dp

  do inp=1,numCells

    CoNum = max( CoNum , su(inp)/Vol(inp) )

    meanCoNum = meanCoNum + su(inp)
    
    suma = suma + Vol(inp)

  enddo

  ! AllReduce using sum, values of suma and MeanCoNum.
  call global_sum(MeanCoNum)
  call global_sum(suma)

  CoNum = 0.5*CoNum*timestep

  ! Find global maximum Courant number in the whole field.
  call global_max(CoNum)


  ! Now use it to calculate mean Courant number
  meanCoNum = 0.5*meanCoNum/suma*timestep

  !// If we keep the value of Courant Number fixed
  if( CoNumFix .and. itime.ne.itimes ) then
      dt = timestep
      timestep = CoNumFixValue/CoNum * timestep

      CoNum = CoNum * timestep/dt
      meanCoNum  = meanCoNum * timestep/dt
  endif

  time = time + timestep

if(myid .eq. 0) then
 write(6,*)
 write(6,'(a,i0,a,es10.3,a,es10.3)') "  Time step no. : ",ITIME," dt : ",timestep," Time = ",time
 write(6,*)

 write(6,'(2(a,es10.3))') "  Courant Number mean: ", meanCoNum," max: ", CoNum
endif


endif
!// ************************************************************************* //
