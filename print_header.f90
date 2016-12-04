!***********************************************************************
!
subroutine print_header
!
!***********************************************************************
!
  use types
  use parameters
  use indexes
  use title_mod

  implicit none 
!
!***********************************************************************
!
  WRITE(6,'(//,45X,A)') '######################################################################'
  WRITE(6,'(50X,A70)')                             TITLE
  WRITE(6,'(45X,A)')    '######################################################################'
  WRITE(6,'(/,50X,A,1PE10.4)') 'FLOW INLET       : FLOW = ',FLOMAS
  WRITE(6,'(50X,A,1PE10.4)')   'FLUID DENSITY    :  DEN = ',DENSIT
  WRITE(6,'(50X,A,1PE10.4)')   'DYNAMIC VISCOSITY:  VIS = ',VISCOS
  WRITE(6,'(50X,A,1PE10.4)')   'CONV. CRITERION  :  SOR = ',SORMAX
  !------------------------------------------------------------------------------------------------
  IF(LCAL(IU))   WRITE(6,'(/,50X,A,F5.2,A,F5.2,A,I3)') 'URF( U )=',URF(IU), '  GDS( U )=',GDS(IU), '  NSW=', NSW(IU)
  IF(LCAL(IV))   WRITE(6,'(50X,A,F5.2,A,F5.2,A,I3)')   'URF( V )=',URF(IV), '  GDS( V )=',GDS(IU), '  NSW=', NSW(IV)
  IF(LCAL(IW))   WRITE(6,'(50X,A,F5.2,A,F5.2,A,I3)')   'URF( W )=',URF(IW), '  GDS( W )=',GDS(IU), '  NSW=', NSW(IW)
  IF(LCAL(IP))   WRITE(6,'(50X,A,F5.2,15X,A,I3)')      'URF( P )=',URF(IP),                       '   NSW=', NSW(IP)
  IF(LCAL(ITE))  WRITE(6,'(50X,A,F5.2,A,F5.2,A,I3)')   'URF( TE)=',URF(ITE),'  GDS( TE)=',GDS(ITE),'  NSW=', NSW(ITE)
  IF(LCAL(IED))  WRITE(6,'(50X,A,F5.2,A,F5.2,A,I3)')   'URF( ED)=',URF(IED),'  GDS( ED)=',GDS(IED),'  NSW=', NSW(IED)
  IF(LCAL(IVIS)) WRITE(6,'(50X,A,F5.2)')               'URF(VIS)=',URF(IVIS)
  IF(LCAL(IEN))  WRITE(6,'(50X,A,F5.2,A,F5.2,A,I3)')   'URF( T )=',URF(IEN),'  GDS( T )=',GDS(IEN),'  NSW=',NSW(IEN)
  IF(LCAL(IVART))WRITE(6,'(50X,A,F5.2,A,F5.2,A,I3)')   'URF(VAR)=',URF(IVART),'  GDS(VAR)=',GDS(IVART),'  NSW=',NSW(IVART)
  IF(LCAL(ICON)) WRITE(6,'(50X,A,F5.2,A,F5.2,A,I3)')   'URF(CON)=',URF(ICON),'  GDS(CON)=',GDS(ICON),'  NSW=',NSW(ICON)
  !------------------------------------------------------------------------------------------------
  WRITE(6,'(a)') '   '
  WRITE(6,'(45X,A)') '================================================================='
  WRITE(6,'(50X,A,E10.4)') 'TIME STEP= ',TIMESTEP
  WRITE(6,'(45X,A)') '-----------------------------------------------------------------'
  IF(LBUOY) THEN
    WRITE(6,'(50X,A)') 'BUOYANCY ACTIVATED:  - YES-    '
      IF(BOUSSINESQ.EQ.1) THEN
        WRITE(6,'(50X,A)') 'BOUSSINESQ APROXIMATON: - YES- '
      ELSE
        WRITE(6,'(50X,A)') 'BOUSSINESQ APROXIMATON: - NO - '
      ENDIF
  ELSE
    WRITE(6,'(50X,A)') 'BUOYANCY ACTIVATED:  - NO -    '
  ENDIF
  WRITE(6,'(45X,A)') '-----------------------------------------------------------------'
  WRITE(6,'(50X,A,E10.4,2X,E10.4,2X,E10.4)') 'GRAVITY: (X-Y-Z): ',GRAVX,GRAVY,GRAVZ 
  WRITE(6,'(45X,A)') '-----------------------------------------------------------------'
  WRITE(6,'(50X,A,I8)')'No. Cells = ',numCells
  WRITE(6,'(45X,A)') '================================================================='
  WRITE(6,'(a)') '   '
  WRITE(6,'(a)') '   '

end subroutine
