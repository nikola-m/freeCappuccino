      SUBROUTINE CGSTAB_SIP(FI,IFI)
!
!***********************************************************************
!
!    This routine incorporates the CGSTAB solver for seven-diagonal,
!    non-symmetric coefficient matrices (suitable for convection/
!    diffusion problems). See Sect. 5.3.7 for details. Array index
!    IJK calculated from indices I, J, and K according to Table 3.1.
!
!    Writen by Samir Muzaferija, Institut fuer Schiffbau, Hamburg, 1995.
!
!    Original routine modified to incorporate SIP as a preconditioner
!    by Nikola Mirkov, 28.01.2014. nmirkov@vinca.rs
!
!***********************************************************************
!
      USE TYPES
      USE PARAMETERS
      USE INDEXES
      USE COEF
      USE COEFB
      USE TITLE_MOD

      IMPLICIT NONE
!
!***********************************************************************
!
      INTEGER, INTENT(IN) :: IFI
      REAL(PREC), DIMENSION(NXYZ) :: FI 

!
!     LOCAL VARIABLES
!
      INTEGER :: I, J, K, IJK, NS, L
      REAL(PREC), DIMENSION(NXYZ) :: RESO,PK,UK,ZK,VK
      REAL(PREC) :: RES0, RSM, RESMAX, RESL, P1, P2, P3
      REAL(PREC) :: ALF, BETO, GAM, BET, OM, VRES, VV, UKRESO

!.....MAX NO. OF INNER ITERS
      RESMAX = SOR(IFI) 
!
!.....CALCULATE INITIAL RESIDUAL VECTOR
!
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J

            RES(IJK)=SU(IJK)-(AE(IJK)*FI(IJK+NJ) +AW(IJK)*FI(IJK-NJ) &
                             +AN(IJK)*FI(IJK+1)  +AS(IJK)*FI(IJK-1) &
                             +AT(IJK)*FI(IJK+NIJ)+AB(IJK)*FI(IJK-NIJ)) &
                            -AP(IJK)*FI(IJK)
          END DO
        END DO
      END DO

      DO I=1,NOC
        RES(IJL(I))=RES(IJL(I))-AR(I)*FI(IJR(I))
        RES(IJR(I))=RES(IJR(I))-AL(I)*FI(IJL(I))
      END DO

      RES0=0.0D0
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            RES0=RES0+ABS(RES(IJK))
          END DO
        END DO
      END DO
!
      IF(LTEST) WRITE(66,*) '                   0  SWEEP, RES0 = ',RES0

!.....CALCULATE COEFFICIENTS OF  L  AND  U  MATRICES USING SIP
      DO K=2,NKM
      DO I=2,NIM
      DO J=2,NJM
       IJK=LK(K)+LI(I)+J
          BB(IJK)=AB(IJK)/(1.+ALFA*(BN(IJK-NIJ)+BE(IJK-NIJ)))
          BW(IJK)=AW(IJK)/(1.+ALFA*(BN(IJK-NJ)+BT(IJK-NJ)))
          BS(IJK)=AS(IJK)/(1.+ALFA*(BE(IJK-1)+BT(IJK-1)))
          P1=ALFA*(BB(IJK)*BN(IJK-NIJ)+BW(IJK)*BN(IJK-NJ))
          P2=ALFA*(BB(IJK)*BE(IJK-NIJ)+BS(IJK)*BE(IJK-1))
          P3=ALFA*(BW(IJK)*BT(IJK-NJ)+BS(IJK)*BT(IJK-1))
          BP(IJK)=1./(AP(IJK)+P1+P2+P3 &
            -BB(IJK)*BT(IJK-NIJ) &
            -BW(IJK)*BE(IJK-NJ) &
            -BS(IJK)*BN(IJK-1)+SMALL)
           BN(IJK)=(AN(IJK)-P1)*BP(IJK)
           BE(IJK)=(AE(IJK)-P2)*BP(IJK)
           BT(IJK)=(AT(IJK)-P3)*BP(IJK)
          END DO
        END DO
      END DO
!
!.....INITIALIZE WORKING ARRAYS AND CONSTANTS
!
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            RESO(IJK)=RES(IJK)
            PK(IJK)=0.0D0
            UK(IJK)=0.0D0
            ZK(IJK)=0.0D0
            VK(IJK)=0.0D0
          END DO
        END DO
      END DO
      ALF=1.0D0
      BETO=1.0D0
      GAM=1.0D0
!
!.....START INNER ITERATIONS
!
      NS=NSW(IFI)
      DO L=1,NS
!
!..... CALCULATE BETA AND OMEGA
!
      BET=0.0D0
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            BET=BET+RES(IJK)*RESO(IJK)
          END DO
        END DO
      END DO
      OM=BET*GAM/(ALF*BETO+SMALL)
      BETO=BET
!
!..... CALCULATE PK
!
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            PK(IJK)=RES(IJK)+OM*(PK(IJK)-ALF*UK(IJK))
          END DO
        END DO
      END DO
!
!.....SOLVE (M ZK = PK) - FORWARD SUBSTITUTION
!
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=(PK(IJK)-BB(IJK)*ZK(IJK-NIJ)-BW(IJK)*ZK(IJK-NJ)- &
            BS(IJK)*ZK(IJK-1))*BP(IJK)
          END DO
        END DO
      END DO

!
!..... BACKWARD ELIMINATION
!
      DO K=NKM,2,-1
        DO I=NIM,2,-1
          DO J=NJM,2,-1
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=ZK(IJK)-BN(IJK)*ZK(IJK+1)-BE(IJK)*ZK(IJK+NJ)- &
                    BT(IJK)*ZK(IJK+NIJ)
          END DO
        END DO
      END DO
!
!.....CALCULATE UK = A.PK
!
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            UK(IJK)=AP(IJK)*ZK(IJK)-(AE(IJK)*ZK(IJK+NJ) +AW(IJK)*ZK(IJK-NJ) &
                                    +AN(IJK)*ZK(IJK+1)  +AS(IJK)*ZK(IJK-1) &
                                    +AT(IJK)*ZK(IJK+NIJ)+AB(IJK)*ZK(IJK-NIJ))
          END DO
        END DO
      END DO
!
!..... CALCULATE SCALAR PRODUCT UK.RESO AND GAMMA
!
      UKRESO=0.0D0
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            UKRESO=UKRESO+UK(IJK)*RESO(IJK)
          END DO
        END DO
      END DO
      GAM=BET/UKRESO
!
!.....UPDATE (FI) AND CALCULATE W (OVERWRITE RES - IT IS RES-UPDATE)
!
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            FI(IJK)=FI(IJK)+GAM*ZK(IJK)   
            RES(IJK)=RES(IJK)-GAM*UK(IJK) !W
          END DO
        END DO
      END DO
!
!.....SOLVE (M Y = W); Y OVERWRITES ZK; FORWARD SUBSTITUTION
!
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=(RES(IJK)-BB(IJK)*ZK(IJK-NIJ)-BW(IJK)*ZK(IJK-NJ)- &
            BS(IJK)*ZK(IJK-1))*BP(IJK)
          END DO
        END DO
      END DO
!
!.....BACKWARD SUBSTITUTION
!
      DO K=NKM,2,-1
        DO I=NIM,2,-1
          DO J=NJM,2,-1
            IJK=LK(K)+LI(I)+J
            ZK(IJK)=ZK(IJK)-BN(IJK)*ZK(IJK+1)-BE(IJK)*ZK(IJK+NJ)- &
                    BT(IJK)*ZK(IJK+NIJ)
          END DO
        END DO
      END DO
!
!.....CALCULATE V = A.Y (VK = A.ZK)
!
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            VK(IJK)=AP(IJK)*ZK(IJK)-(AE(IJK)*ZK(IJK+NJ) +AW(IJK)*ZK(IJK-NJ) &
                                    +AN(IJK)*ZK(IJK+1)  +AS(IJK)*ZK(IJK-1) &
                                    +AT(IJK)*ZK(IJK+NIJ)+AB(IJK)*ZK(IJK-NIJ))
          END DO
        END DO
      END DO
!
!..... CALCULATE ALPHA (ALF)
!
      VRES=0.0D0
      VV=0.0D0
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            VRES=VRES+VK(IJK)*RES(IJK)
            VV=VV+VK(IJK)*VK(IJK)
          END DO
        END DO
      END DO

      ALF=VRES/(VV+SMALL)
!
!.....UPDATE VARIABLE (FI) AND RESIDUAL (RES) VECTORS
!
      RESL=0.0D0
      DO K=2,NKM
        DO I=2,NIM
          DO J=2,NJM
            IJK=LK(K)+LI(I)+J
            FI(IJK)=FI(IJK)+ALF*ZK(IJK)
            RES(IJK)=RES(IJK)-ALF*VK(IJK)
            RESL=RESL+ABS(RES(IJK))
          END DO
        END DO
      END DO
!
!.....CHECK CONVERGENCE
!
      IF(L.EQ.1) RESOR(IFI)=RES0
      RSM=RESL/(RESOR(IFI)+SMALL)
      IF(LTEST) WRITE(66,'(19x,3a,I4,a,1PE10.3,a,1PE10.3)') ' FI=',CHVAR(IFI),' SWEEP = ',L,' RESL = ',RESL,' RSM = ',RSM
      IF(RSM.LT.RESMAX) EXIT
!
!.....END OF ITERATION LOOP
!
      END DO

!.....Write linear solver report:
      write(66,'(3a,i1,a,i1,2a,1PE10.3,a,1PE10.3,a,I0)') &
      'BiCGStab(SIP):  Solving for ',trim(chvarSolver(IFI)),'(icorr=',icorr,',npcor=',ipcorr,')', &
      ', Initial residual = ',RES0,', Final residual = ',RESL,', No Iterations ',L
!
      END SUBROUTINE

