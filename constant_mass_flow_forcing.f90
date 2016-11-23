         !# Correct driving force for a constant mass flow rate.......................

         ! # Extract the velocity in the flow direction
         ! magUbarStar = ( flowDirection & U ).weightedAverage( mesh.V() )
         magUbarStar = volumeWeightedAverage(U)

         ! # Calculate the pressure gradient increment needed to
         ! # adjust the average flow-rate to the correct value
         ! gragPplus = ( magUbar - magUbarStar ) / rUA.weightedAverage( mesh.V() )
         rUAw = volumeWeightedAverage(APU)
         gragPplus = ( magUbar - magUbarStar ) / rUAw

         ! #  Correction
         ! U.ext_assign( U + flowDirection * rUA * gragPplus )
         flowDirection = 1.0_dp
         U =  U + flowDirection * APU * gragPplus

         ! # Pressure gradient force that will drive the flow - we use it in calcuvw.
         gradPcmf  = gradPcmf + gragPplus

         write(66,'(2(a,es13.6))') "Uncorrected Ubar = ",magUbarStar," pressure gradient = ",gradPcmf
         !............................................................................