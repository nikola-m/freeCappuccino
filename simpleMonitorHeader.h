
      IF(ITIME.EQ.ITIMES) THEN
      WRITE(66,"(/,30X,60(1H=),/,30X,A,I4,A,2X,A,I3,A,I3,A,I3,A,I8,A,/30X,60(1H=))")&
      '* GRID NO ***',KGRID,' ***','N0 CV = ',NI-2,' X ',NJ-2,' X ',NK-2,' = ',(NI-2)*(NJ-2)*(NK-2),' CV *'
      WRITE(66,"(/,1X,A,2X,A,I3,A,I3,A,I3,A)")&
      '                         I---------------ABSOLUTE RESIDUAL SOURCE SUMS---------------I',&
      '       I---------FIELD VALUES AT MONITORING LOCATION (',IIM,',',JJM,',',KKM,')---------I'
      WRITE(66,"(/,1X,A,5X,A,5X,A,5X,A,5X,A,5X,A,5X,A,5X,A,5X,A,5X,A,8X,A,9X,A,9X,A,9X,A,9X,A,7X,A,7X,A,7X,A,7X,A/)")& 
      'ITER','UMOM','VMOM','WMOM','MASS','T EN','DISS','ENER','VART',' CON','U','V','W','P','TE','DISS','T','VART','CON'
      END IF
