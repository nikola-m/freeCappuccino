!***********************************************************************
!
      SUBROUTINE SIPSOL(FI,IFI)
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
      INTEGER :: I, J, K, IJK, NS, N, LKI
      REAL(PREC) :: RSM, RES1, P1, P2, P3

!.....CALCULATE COEFFICIENTS OF  L  AND  U  MATRICES
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
      ENDDO
      ENDDO
      ENDDO
!
!.....ITERATION LOOP
      NS=NSW(IFI)
      DO N=1,NS

!.....CALCULATE RESIDUALS AND AUXILLIARY VECTOR
      RES1=0.0d0

      DO K=2,NKM
      DO I=2,NIM
      LKI=LK(K)+LI(I)
      DO IJK=LKI+2,LKI+NJM
      RES(IJK)=SU(IJK)-(AE(IJK)*FI(IJK+NJ) +AW(IJK)*FI(IJK-NJ) &
                       +AN(IJK)*FI(IJK+1)  +AS(IJK)*FI(IJK-1) &
                       +AT(IJK)*FI(IJK+NIJ)+AB(IJK)*FI(IJK-NIJ)) &
                      -AP(IJK)*FI(IJK)
      ENDDO
      ENDDO
      ENDDO

      DO I=1,NOC
        RES(IJL(I))=RES(IJL(I))-AR(I)*FI(IJR(I))
        RES(IJR(I))=RES(IJR(I))-AL(I)*FI(IJL(I))
      END DO

      DO K=2,NKM
      DO I=2,NIM
      LKI=LK(K)+LI(I)
      DO IJK=LKI+2,LKI+NJM

      RES1=RES1+ABS(RES(IJK))

      RES(IJK)=(RES(IJK)-BB(IJK)*RES(IJK-NIJ)-BW(IJK)*RES(IJK-NJ)- &
            BS(IJK)*RES(IJK-1))*BP(IJK)
      ENDDO
      ENDDO
      ENDDO
      
      IF(N.EQ.1) RESOR(IFI)=RES1
      RSM=RES1/(RESOR(IFI)+SMALL)
!.....CALCULATE INCREMENT AND UPDATE VARIABLES
      DO K=NKM,2,-1
      DO I=NIM,2,-1
      LKI=LK(K)+LI(I)
      DO IJK=LKI+NJM,LKI+2,-1
      RES(IJK)=RES(IJK)-BN(IJK)*RES(IJK+1)-BE(IJK)*RES(IJK+NJ)- &
             BT(IJK)*RES(IJK+NIJ)
      FI(IJK)=FI(IJK)+RES(IJK)
      ENDDO
      ENDDO
      ENDDO
      IF(LTEST) WRITE(66,"(20X,'FI=',A3,'  SWEEP=',I4,'  RES=',1PE10.3,' RSM=',1PE10.3)") CHVAR(IFI),N,RES1,RSM
      IF(RSM.LT.SOR(IFI)) EXIT
      ENDDO
!.....Write linear solver report:
      write(66,'(3a,1PE10.3,a,1PE10.3,a,I0)') &
      'SIP:  Solving for ',trim(chvarSolver(IFI)), &
      ', Initial residual = ',RESOR(IFI),', Final residual = ',RES1,', No Iterations ',N

      END SUBROUTINE
