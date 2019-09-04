!-------------------------------------------------------------------------------
!    initialzie constant variable
!-------------------------------------------------------------------------------

SUBROUTINE INIT_VARIABLES
INCLUDE 'INCLUDE_COMMON.FI'

LOGICAL IEXIST
INTEGER I , J , K 
      

DO I = 1 , NDX
    KX(I) = ARF * DBLE(I-1)
ENDDO

DO J = 1 , NDZ/2
    KZ(J) = BET * DBLE(J-1)
ENDDO

DO J = NDZ/2 + 1 , NDZ
        KZ(J) = BET * DBLE(J - 1 - NDZ)
ENDDO

DO K = 1 , NDT
        OMEGA(K) = GAMMA * DBLE(K - 1 - NDT/2)
ENDDO

IF(IMPI_MYID .EQ. 0) THEN
    INQUIRE(FILE='SPECTRA_T/..',EXIST=IEXIST)
    IF(IEXIST .NEQV. .TRUE.) CALL SYSTEM('mkdir SPECTRA_T')
!-------------------------------------------------------------------------------
!    Read the spectrum files
!-------------------------------------------------------------------------------	
    OPEN (21 , FILE='uu39.dat')
        DO I = 1 , NDX
        DO J = 1 , NDZ
        DO K = 1 , NDT
                READ(21,*)EPUU(I,J,K)
        ENDDO
        ENDDO
        ENDDO
    CLOSE(21)

        OPEN (21 , FILE='vv39.dat')
        DO I = 1 , NDX
        DO J = 1 , NDZ
        DO K = 1 , NDT
                READ(21,*)EPVV(I,J,K)
        ENDDO
        ENDDO
        ENDDO
    CLOSE(21)

        OPEN (21 , FILE='ww39.dat')
        DO I = 1 , NDX
        DO J = 1 , NDZ
        DO K = 1 , NDT
                READ(21,*)EPWW(I,J,K)
        ENDDO
        ENDDO
        ENDDO
    CLOSE(21)
        
        OPEN (21 , FILE='uv39.dat')
        DO I = 1 , NDX
        DO J = 1 , NDZ
        DO K = 1 , NDT
                READ(21,*)EPUV(I,J,K)
        ENDDO
        ENDDO
        ENDDO
    CLOSE(21)

        OPEN (21 , FILE='uw39.dat')
        DO I = 1 , NDX
        DO J = 1 , NDZ
        DO K = 1 , NDT
                READ(21,*)EPUW(I,J,K)
        ENDDO
        ENDDO
        ENDDO
    CLOSE(21)

        OPEN (21 , FILE='vw39.dat')
        DO I = 1 , NDX
        DO J = 1 , NDZ
        DO K = 1 , NDT
                READ(21,*)EPVW(I,J,K)
        ENDDO
        ENDDO
        ENDDO
    CLOSE(21)
!-------------------------------------------------------------------------------
!    Read the Chebyshev files
!-------------------------------------------------------------------------------
        OPEN (21 , FILE='cheb_point.dat')
        DO I = 1 , NDY
                READ(21,*)CHEB_POINT(I)
        ENDDO
    CLOSE(21)

        OPEN (21 , FILE='cheb_MAT1.dat')
        DO I = 1 , NDY
        DO J = 1 , NDY
                READ(21,*)CHEB_MAT1(I,J)
        ENDDO
        ENDDO
    CLOSE(21)

        OPEN (21 , FILE='cheb_MAT2.dat')
        DO I = 1 , NDY
        DO J = 1 , NDY
                READ(21,*)CHEB_MAT2(I,J)
        ENDDO
        ENDDO
    CLOSE(21)

!-------------------------------------------------------------------------------
!    Read the Profile files
!-------------------------------------------------------------------------------
        OPEN (21 , FILE='U.dat')
        DO I = 1 , NDY
                READ(21,*)U(I)
        ENDDO
    CLOSE(21)

        OPEN (21 , FILE='DU.dat')
        DO I = 1 , NDY
                READ(21,*)DU(I)
        ENDDO
    CLOSE(21)

        OPEN (21 , FILE='NU_T.dat')
        DO I = 1 , NDY
                READ(21,*)NU_T(I)
        ENDDO
    CLOSE(21)

        OPEN (21 , FILE='DNU_T.dat')
        DO I = 1 , NDY
                READ(21,*)DNU_T(I)
        ENDDO
    CLOSE(21)

ENDIF

CALL MPI_BCAST(EPUU         , NDX*NDZ*NDT , MPI_DOUBLE_COMPLEX   ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)
CALL MPI_BCAST(EPVV         , NDX*NDZ*NDT , MPI_DOUBLE_COMPLEX   ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)
CALL MPI_BCAST(EPWW         , NDX*NDZ*NDT , MPI_DOUBLE_COMPLEX   ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)
CALL MPI_BCAST(EPUV         , NDX*NDZ*NDT , MPI_DOUBLE_COMPLEX   ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)
CALL MPI_BCAST(EPUW         , NDX*NDZ*NDT , MPI_DOUBLE_COMPLEX   ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)
CALL MPI_BCAST(EPVW         , NDX*NDZ*NDT , MPI_DOUBLE_COMPLEX   ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)

CALL MPI_BCAST(CHEB_POINT   , NDY         , MPI_DOUBLE_PRECISION ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)
CALL MPI_BCAST(CHEB_MAT1    , NDY*NDY     , MPI_DOUBLE_PRECISION ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)
CALL MPI_BCAST(CHEB_MAT2    , NDY*NDY     , MPI_DOUBLE_PRECISION ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)

CALL MPI_BCAST(U            , NDY         , MPI_DOUBLE_PRECISION ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)
CALL MPI_BCAST(DU           , NDY         , MPI_DOUBLE_PRECISION ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)
CALL MPI_BCAST(NU_T         , NDY         , MPI_DOUBLE_PRECISION ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)
CALL MPI_BCAST(DNU_T        , NDY         , MPI_DOUBLE_PRECISION ,  0 , MPI_COMM_WORLD , iMPI_ErrorInfo)


H         = 0.
RQ        = 0.
RY        = 0.
Us        = 0.
Vt        = 0.
S         = 0.
SIGMA     = 0.
INV_SIGMA = 0.
V1        = 0.
TRAN_V1   = 0.
TRAN_U    = 0.
TRAN_RQ   = 0.
SYY       = 0.
SQQ       = 0.
ET        = 0.
ETT       = 0.
E11       = 0.
TEMP1     = 0.
TEMP2     = 0.

U_ES1     = 0.
U_ES2     = 0.
U_ES3     = 0.
U_ES4     = 0.


END
