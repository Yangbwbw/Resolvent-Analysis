!----------------------------------------------------------------------------------------
! Subroutine for writing files of spatio-temporal spectrum
!----------------------------------------------------------------------------------------

SUBROUTINE WRITEFILE
INCLUDE 'INCLUDE_COMMON.FI'
INTEGER I , J , K , L
CHARACTER*3 YNAME
INTEGER FNFULL , FNF
INTEGER FNI(3)
CHARACTER FNC(3)

FNF = INT(IMPI_MYID)
FNFULL = FNF
        DO I = 3 , 1 , -1
                FNI(I) = MOD(FNFULL , 10)
                FNFULL = (FNFULL-FNI(I))/10
                WRITE(FNC(I) , '(I1)') FNI(I)
        ENDDO
        YNAME = FNC(1)//FNC(2)//FNC(3)
        OPEN(21 , FILE = 'SPECTRA_T/'//'ET'//YNAME//'.DAT')
        DO J = 1 , 3
        DO K = 1 , NET
                WRITE(21 , *)ET(J , K)
        ENDDO
        ENDDO
        CLOSE(21)
		
		OPEN(21 , FILE = 'SPECTRA_T/'//'U_ES1'//YNAME//'.DAT')
        DO I = 1 , 3
        DO K = 1 , NDTNP
		DO J = 1 , NDX
                WRITE(21 , *)U_ES1(J , K , I)
        ENDDO
        ENDDO
		ENDDO
        CLOSE(21)
		
		OPEN(21 , FILE = 'SPECTRA_T/'//'U_ES2'//YNAME//'.DAT')
        DO I = 1 , 3
        DO K = 1 , NDTNP
		DO J = 1 , NDX
                WRITE(21 , *)U_ES2(J , K , I)
        ENDDO
        ENDDO
		ENDDO
        CLOSE(21)
		
		OPEN(21 , FILE = 'SPECTRA_T/'//'U_ES3'//YNAME//'.DAT')
        DO I = 1 , 3
        DO K = 1 , NDTNP
		DO J = 1 , NDX
                WRITE(21 , *)U_ES3(J , K , I)
        ENDDO
        ENDDO
		ENDDO
        CLOSE(21)
		
		OPEN(21 , FILE = 'SPECTRA_T/'//'U_ES4'//YNAME//'.DAT')
        DO I = 1 , 3
        DO K = 1 , NDTNP
		DO J = 1 , NDX
                WRITE(21 , *)U_ES4(J , K , I)
        ENDDO
        ENDDO
		ENDDO
        CLOSE(21)

END
