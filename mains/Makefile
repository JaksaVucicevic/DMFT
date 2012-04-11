# C++ compiler
Cpp =icpc #/opt/intel/bin/icpc

main =main_tmt

#mpiCC
# SIAM path
SP = ../source
# executable path
RP = run

# -----------  SERIAL ------------------------------------#

#mpiCC = icpc
#FLAGS =  -static-intel

# ----------- THREADED -----------------------------------# 

mpiCC = icpc
FLAGS =  -D_OMP -openmp #-static-intel

# -----------  HYBRID ------------------------------------#

#mpiCC = mpiCC
#FLAGS =  -D_MPI -D_OMP -openmp -static-intel #-fast

#---------------------------------------------------------#

LIBS =# use this if needed 

all : $(main).o $(SP)/TMT.o $(SP)/CHM.o $(SP)/Loop.o $(SP)/SIAM.o $(SP)/Result.o $(SP)/GRID.o $(SP)/Input.o $(SP)/Broyden.o $(SP)/Broyden.h $(SP)/Mixer.h $(SP)/routines.o $(SP)/nrutil.o
	$(mpiCC) $(FLAGS) -o $(RP)/$(main) $(LIBS) $(main).o $(SP)/TMT.o $(SP)/CHM.o $(SP)/Loop.o $(SP)/SIAM.o $(SP)/Result.o $(SP)/GRID.o $(SP)/Input.o $(SP)/Broyden.o $(SP)/routines.o $(SP)/nrutil.o

# main program
$(main).o : $(main).cpp $(SP)/TMT.h $(SP)/CHM.h $(SP)/SIAM.h $(SP)/Result.h $(SP)/GRID.h
	$(mpiCC) $(FLAGS) -c -o $@ $(main).cpp

# TMT
$(SP)/TMT.o : $(SP)/TMT.cpp $(SP)/TMT.h $(SP)/CHM.h $(SP)/Loop.h $(SP)/SIAM.h $(SP)/Result.h $(SP)/GRID.h $(SP)/Input.h $(SP)/routines.h
	$(mpiCC) $(FLAGS) -c -o $@ $(SP)/TMT.cpp

# CHM
$(SP)/CHM.o : $(SP)/CHM.cpp $(SP)/CHM.h $(SP)/Loop.h $(SP)/SIAM.h $(SP)/Result.h $(SP)/GRID.h $(SP)/Input.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/CHM.cpp

# Loop (base class for CHM and TMT)
$(SP)/Loop.o : $(SP)/Loop.h $(SP)/Loop.cpp $(SP)/Result.h $(SP)/GRID.h $(SP)/Input.h $(SP)/Mixer.h $(SP)/Broyden.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Loop.cpp

# SIAM
$(SP)/SIAM.o : $(SP)/SIAM.cpp $(SP)/SIAM.h $(SP)/Broyden.h $(SP)/Result.h $(SP)/GRID.h $(SP)/Input.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/SIAM.cpp

# Result
$(SP)/Result.o : $(SP)/Result.cpp $(SP)/Result.h $(SP)/GRID.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Result.cpp

# Grid utility for initializing omega grids and provides all grid dependent routines
$(SP)/GRID.o : $(SP)/GRID.cpp $(SP)/GRID.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/GRID.cpp

# Input class used for reading files with parameters
$(SP)/Input.o : $(SP)/Input.cpp $(SP)/Input.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Input.cpp

# Broyden solver for systems of equations
$(SP)/Broyden.o : $(SP)/Broyden.h $(SP)/Broyden.cpp
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Broyden.cpp

# contains some constants and useful numerical routines
$(SP)/routines.o : $(SP)/routines.cpp $(SP)/routines.h 
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/routines.cpp

# numerical routines from NumRec
$(SP)/nrutil.o : $(SP)/nrutil.h $(SP)/nrutil.c
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/nrutil.c

# clean all object and exec files
clean :
	rm -f $(RP)/$(main) *.o $(SP)/*.o
