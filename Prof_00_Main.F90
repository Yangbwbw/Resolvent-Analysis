Program Resolvent_Analysis
INCLUDE 'INCLUDE_COMMON.FI'

INTEGER I , J , K , L

CALL INIT_PARALLEL
CALL INIT_VARIABLES

DO II = 1 , NDX
DO KK = 1 , NDZ
DO TT = TT_S , TT_E

	IF ((II .ne. 1).or.(KK .ne. 1)) THEN

        KXX = KX(II)
        KZZ = KZ(KK)
        OMEGAA = OMEGA(TT)
        K2 = KXX*KXX + KZZ*KZZ

        CALL RESET_VARIABLES
!-------------------------------------------------------------------------------
!    Calculate the resolvent matrix
!-------------------------------------------------------------------------------
                DO J = 1 , NDY
                        H(4*J-3 , 4*J-3) = -IC * OMEGAA + IC * KXX * U(J) + K2 * NU_T(J) / RE
                        H(4*J   , 4*J-3) = IC * KXX
                        H(4*J-3 , 4*J-2) = DU(J) - IC * KXX / RE * DNU_T(J)
                        H(4*J-2 , 4*J-2) = -IC * OMEGAA + IC * KXX * U(J) + K2 * NU_T(J) / RE
                        H(4*J-1 , 4*J-2) = -IC * KZZ / RE * DNU_T(J)
                        H(4*J-1 , 4*J-1) = -IC * OMEGAA + IC * KXX * U(J) + K2 * NU_T(J) / RE
                        H(4*J   , 4*J-1) = IC * KZZ
                        H(4*J-3 , 4*J  ) = IC * KXX
                        H(4*J-1 , 4*J  ) = IC * KZZ
                ENDDO

                DO J = 1 , NDY
                DO K = 1 , NDY
                        H(4*J-3 , 4*K-3) = H(4*J-3 , 4*K-3) - 1.0 / RE * DNU_T(J) * CHEB_MAT1(J , K) - NU_T(J) / RE * CHEB_MAT2(J , K)
                        H(4*J-2 , 4*K-2) = H(4*J-2 , 4*K-2) - 2.0 / RE * DNU_T(J) * CHEB_MAT1(J , K) - NU_T(J) / RE * CHEB_MAT2(J , K)
                        H(4*J-1 , 4*K-1) = H(4*J-1 , 4*K-1) - 1.0 / RE * DNU_T(J) * CHEB_MAT1(J , K) - NU_T(J) / RE * CHEB_MAT2(J , K)
                        H(4*J-2 , 4*K  ) = H(4*J-2 , 4*K  ) + CHEB_MAT1(J , K)
                        H(4*J   , 4*K-2) = H(4*J   , 4*K-2) + CHEB_MAT1(J , K)
                ENDDO
                ENDDO

!-------------------------------------------------------------------------------
!    Velocity boundary condition
!-------------------------------------------------------------------------------
                ! H(4    , :) = H(2    , :)
				! H(NH   , :) = H(NH-2 , :)
				
				H(1    , :) = 0.
                H(2    , :) = 0.
                H(3    , :) = 0.
                H(NH-1 , :) = 0.
                H(NH-2 , :) = 0.
                H(NH-3 , :) = 0.
                
                H(1 , 1) = 1.0
                H(2 , 2) = 1.0
                H(3 , 3) = 1.0
                H(NH-1 , NH-1) = 1.0
                H(NH-2 , NH-2) = 1.0
                H(NH-3 , NH-3) = 1.0
				
!-------------------------------------------------------------------------------
!    Calculate the resolvent
!-------------------------------------------------------------------------------
                CALL INV_MAT(NH , H , R)

!-------------------------------------------------------------------------------
!    Delete the columns of continuity equations 
!-------------------------------------------------------------------------------
                DO J = 2 , NDY-1
                        RQ(: , 3*(J-2)+1) = R(: , 4*(J-1)+1)
                        RQ(: , 3*(J-2)+2) = R(: , 4*(J-1)+2)
                        RQ(: , 3*(J-2)+3) = R(: , 4*(J-1)+3)
                ENDDO

                        RY(1 , :)         = RQ(4*PG+1 , :)
                        RY(2 , :)         = RQ(4*PG+2 , :)
                        RY(3 , :)         = RQ(4*PG+3 , :)
!-------------------------------------------------------------------------------
!    Transpose of RQ
!-------------------------------------------------------------------------------
                DO J = 1 , NH
                DO K = 1 , MRY
                        TRAN_RQ(K , J) = DCONJG(RQ(J , K))
                ENDDO
                ENDDO

!-------------------------------------------------------------------------------
!    Singular Value Decomposition(SVD)
!-------------------------------------------------------------------------------

                CALL SVD(NRY , MRY , RY , S , Us , Vt)

!-------------------------------------------------------------------------------
!    Transpose of Singular Vector
!-------------------------------------------------------------------------------
                DO J = 1 , NRY
                DO K = 1 , NRY
                        TRAN_U(K , J) = DCONJG(Us(J , K))
                ENDDO
                ENDDO

                DO J = 1 , NRY
                        SIGMA(J , J) = S(J)
                        INV_SIGMA(J , J) = 1.0/S(J)
                ENDDO
                
                DO J = 1 , NRY
                        TRAN_V1(J , :) = Vt(J , :)
                ENDDO

                DO J = 1 , NRY
                DO K = 1 , MRY
                        V1(K , J) = DCONJG(TRAN_V1(J , K))
                ENDDO
                ENDDO

!-------------------------------------------------------------------------------
!    Known Spectra
!-------------------------------------------------------------------------------

                SYY(1 , 1) = EPUU(II , KK , TT)
                SYY(2 , 2) = EPVV(II , KK , TT)
                SYY(3 , 3) = EPWW(II , KK , TT)
                SYY(1 , 2) = EPUV(II , KK , TT)
                SYY(1 , 3) = EPUW(II , KK , TT)
                SYY(2 , 3) = EPVW(II , KK , TT)
                SYY(2 , 1) = DCONJG(EPUV(II , KK , TT))
                SYY(3 , 1) = DCONJG(EPUW(II , KK , TT))
                SYY(3 , 2) = DCONJG(EPVW(II , KK , TT))


                E11 = MATMUL(MATMUL(MATMUL(MATMUL(INV_SIGMA , TRAN_U) , SYY) , Us) , INV_SIGMA)
                SQQ = MATMUL(MATMUL(MATMUL(MATMUL(RQ , V1) , E11) , TRAN_V1) , TRAN_RQ)

                ! DO L = 1 , NET
                        ! ET(1 , L) = ET(1 , L) + DBLE(SQQ(4*L+1 , 4*L+1))
                        ! ET(2 , L) = ET(2 , L) + DBLE(SQQ(4*L+2 , 4*L+2))
                        ! ET(3 , L) = ET(3 , L) + DBLE(SQQ(4*L+3 , 4*L+3))
                ! ENDDO
				IF (II .eq. 1) THEN
					DO L = 1 , NET
							ET(1 , L) = ET(1 , L) + DBLE(SQQ(4*L+1 , 4*L+1))*BET*ARF*GAMMA
							ET(2 , L) = ET(2 , L) + DBLE(SQQ(4*L+2 , 4*L+2))*BET*ARF*GAMMA
							ET(3 , L) = ET(3 , L) + DBLE(SQQ(4*L+3 , 4*L+3))*BET*ARF*GAMMA
					ENDDO
				ELSE
					DO L = 1 , NET
							ET(1 , L) = ET(1 , L) + 2.0*DBLE(SQQ(4*L+1 , 4*L+1))*BET*ARF*GAMMA
							ET(2 , L) = ET(2 , L) + 2.0*DBLE(SQQ(4*L+2 , 4*L+2))*BET*ARF*GAMMA
							ET(3 , L) = ET(3 , L) + 2.0*DBLE(SQQ(4*L+3 , 4*L+3))*BET*ARF*GAMMA
					ENDDO
				ENDIF
				
				DO L = 1 , 3
						U_ES1(II , TT - TT_S +1 , L) = U_ES1(II , TT - TT_S +1 , L) + DBLE(SQQ(4*PE1 + L , 4*PE1 + L))*BET
						U_ES2(II , TT - TT_S +1 , L) = U_ES2(II , TT - TT_S +1 , L) + DBLE(SQQ(4*PE2 + L , 4*PE2 + L))*BET
						U_ES3(II , TT - TT_S +1 , L) = U_ES3(II , TT - TT_S +1 , L) + DBLE(SQQ(4*PE3 + L , 4*PE3 + L))*BET
						U_ES4(II , TT - TT_S +1 , L) = U_ES4(II , TT - TT_S +1 , L) + DBLE(SQQ(4*PE4 + L , 4*PE4 + L))*BET
				ENDDO
	ENDIF

ENDDO
ENDDO
ENDDO

CALL WRITEFILE

CALL FINISH_PARALLEL

END



!--------------------------------------------------------------------------
!     subroutine for calculate the inversion matrix
!--------------------------------------------------------------------------
SUBROUTINE INV_MAT(N , A , INV_A)
INTEGER :: N , LDA , IPIV(N) , INFO , LWORK
COMPLEX*16 :: A(N , N) , INV_A(N , N) , WORK(64 * N)

        INV_A = A
        LWORK = 64 * N
        LDA = N

        CALL zgetrf(N, N, INV_A, LDA, IPIV, INFO)
        IF(INFO .NE. 0) WRITE(0,*) 'Error occured in zgetrf!'
        CALL zgetri(N, INV_A, LDA, IPIV, WORK, LWORK, INFO)
        IF(INFO .NE. 0) WRITE(0,*) 'Error occured in zgetri!'

END

!--------------------------------------------------------------------------
!     subroutine for calculate the Singular Value Decomposition(SVD)
!--------------------------------------------------------------------------
SUBROUTINE SVD(M , N , A , S , Us, Vt)
CHARACTER*1 :: JOBZ
INTEGER :: M , N , LDA , LDU , LDVT , LWORK , IWORK(16 * M) , INFO
REAL*8 :: RWORK(3 * M * N) , S(M)
COMPLEX*16 :: WORK(128 * N) , A(M , N) , Us(M , M) , Vt(N , N)

        JOBZ = 'A'
        LWORK = 128 * N
        LDA = M
        LDU = M
        LDVT = N

        CALL zgesdd(JOBZ , M , N , A , LDA , S , Us , LDU , Vt ,LDVT , WORK , LWORK , RWORK , IWORK , INFO)
        IF(INFO .NE. 0) WRITE(0,*) 'Error occured in zgesdd!'

END
