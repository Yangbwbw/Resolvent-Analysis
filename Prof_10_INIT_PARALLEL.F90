SUBROUTINE INIT_PARALLEL
!-------------------------------------------------------------------------------
!   SUBROUTINE FOR INITIALIZING THE PARALLEL PARAMETERS
!-------------------------------------------------------------------------------
INCLUDE 'INCLUDE_COMMON.FI'

! execute statements
CALL MPI_INIT( iMPI_ErrorInfo )
CALL MPI_COMM_RANK(MPI_COMM_WORLD,IMPI_MYID,iMPI_ErrorInfo)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NUMPROCS,iMPI_ErrorInfo)
!
IF(NUMPROCS.NE.NP) THEN
   IF(IMPI_MYID.EQ.0) WRITE(*,*) 'running process not equal to the process required'
   CALL FINISH_PARALLEL
ENDIF

TT_S = IMPI_MYID * NDTNP + 1
TT_E = TT_S - 1 + NDTNP
END
