### user configurable options #####

FC          = mpiifort
FLINKER     = mpiifort
OPTINCLUDE  = -I/THL7/home/hegw1/Yangbowen/opt/lapack/include
OPTLINK     = -L/THL7/home/hegw1/Yangbowen/opt/lapack/lib -llapack -lrefblas -lz -lm
#### End User cenfigurable options ###

FFLAGS =  $(OPTFLAGS)
FLIBS =
EXECS =  Resolvent
OBJECT =  Prof_00_Main.o\
          Prof_10_INIT_PARALLEL.o\
          Prof_12_INIT_VARIABLES.o\
	  Prof_13_RESET_VARIABLES.o\
          Prof_30_WRITEFILE.o\
          Prof_99_FINALIZE.o
default: $(EXECS)

$(EXECS): $(OBJECT) 
	$(FLINKER) $(OBJECT) $(OPTLINK) -o $(EXECS)

%.o: %.F90
	$(FC) $(OPTINCLUDE) -c $< -o $@

clean:
	/bin/rm -f *.o $(EXECS)
