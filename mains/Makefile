# C++ compiler
Cpp = icpc #/opt/intel/Compiler/11.1/046/bin/intel64/icpc #/opt/intel/bin/icpc

main =main_chm

#mpiCC
# source path
SP = ../source
#object files path
OP = ../obj
# executable path
RP = ../bin

# -----------  SERIAL ------------------------------------#

#mpiCC = icpc #/opt/intel/Compiler/11.1/046/bin/intel64/icpc
#FLAGS =#-static-intel

# ----------- THREADED -----------------------------------# 

mpiCC = icpc #/opt/intel/Compiler/11.1/046/bin/intel64/icpc
FLAGS =   -openmp -D_OMP -static

# -----------  HYBRID ------------------------------------#

#mpiCC = /opt/openmpi-1.6.3/bin/mpicxx
#FLAGS =  -D_MPI -D_OMP -openmp -static-intel #-fast

#---------------------------------------------------------#

LIBS =# use this if needed 

all : $(OP)/$(main).o $(OP)/TMT.o $(OP)/CHM.o $(OP)/CHM2.o $(OP)/Loop.o $(OP)/SIAM.o $(OP)/SIAM2.o $(OP)/LambdaCalculator.o $(OP)/Result.o $(OP)/GRID.o $(OP)/Input.o $(OP)/Broyden.o $(OP)/IBZ.o $(SP)/Broyden.h $(SP)/Mixer.h $(OP)/routines.o $(OP)/nrutil.o
	$(mpiCC) $(FLAGS) -o $(RP)/$(main) $(LIBS) $(OP)/$(main).o  $(OP)/TMT.o $(OP)/CHM.o $(OP)/CHM2.o $(OP)/Loop.o $(OP)/SIAM.o $(OP)/SIAM2.o $(OP)/LambdaCalculator.o $(OP)/Result.o $(OP)/GRID.o $(OP)/Input.o $(OP)/Broyden.o $(OP)/IBZ.o $(OP)/routines.o $(OP)/nrutil.o

# main program
$(OP)/$(main).o : $(main).cpp $(SP)/TMT.h $(SP)/CHM.h $(SP)/CHM2.h $(SP)/SIAM.h $(SP)/SIAM2.h $(SP)/Result.h $(SP)/GRID.h
	$(mpiCC) $(FLAGS) -c -o $@ $(main).cpp

# TMT
$(OP)/TMT.o : $(SP)/TMT.cpp $(SP)/TMT.h $(SP)/CHM.h $(SP)/Loop.h $(SP)/SIAM.h $(SP)/Result.h $(SP)/GRID.h $(SP)/Input.h $(SP)/routines.h
	$(mpiCC) $(FLAGS) -c -o $@ $(SP)/TMT.cpp

# CHM
$(OP)/CHM.o : $(SP)/CHM.cpp $(SP)/CHM.h $(SP)/Loop.h $(SP)/SIAM.h $(SP)/Result.h $(SP)/GRID.h $(SP)/Input.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/CHM.cpp

# CHM2
$(OP)/CHM2.o : $(SP)/CHM2.cpp $(SP)/CHM2.h $(SP)/Loop.h $(SP)/SIAM2.h $(SP)/Result.h $(SP)/GRID.h $(SP)/Input.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/CHM2.cpp

# Loop (base class for CHM and TMT)
$(OP)/Loop.o : $(SP)/Loop.h $(SP)/Loop.cpp $(SP)/Result.h $(SP)/GRID.h $(SP)/Input.h $(SP)/Mixer.h $(SP)/Broyden.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Loop.cpp

# SIAM
$(OP)/SIAM.o : $(SP)/SIAM.cpp $(SP)/SIAM.h $(SP)/Broyden.h $(SP)/Result.h $(SP)/GRID.h $(SP)/Input.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/SIAM.cpp

# SIAM2
$(OP)/SIAM2.o : $(SP)/SIAM2.cpp $(SP)/SIAM2.h $(SP)/Broyden.h $(SP)/Result.h $(SP)/GRID.h $(SP)/Input.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/SIAM2.cpp

# LambdaCalculator
$(OP)/LambdaCalculator.o : $(SP)/LambdaCalculator.cpp $(SP)/LambdaCalculator.h $(SP)/Input.h $(SP)/routines.h 
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/LambdaCalculator.cpp


# Result
$(OP)/Result.o : $(SP)/Result.cpp $(SP)/Result.h $(SP)/GRID.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Result.cpp

# Grid utility for initializing omega grids and provides all grid dependent routines
$(OP)/GRID.o : $(SP)/GRID.cpp $(SP)/GRID.h $(SP)/routines.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/GRID.cpp

# Input class used for reading files with parameters
$(OP)/Input.o : $(SP)/Input.cpp $(SP)/Input.h
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Input.cpp

# Broyden solver for systems of equations
$(OP)/Broyden.o : $(SP)/Broyden.h $(SP)/Broyden.cpp
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/Broyden.cpp

# IBZ
$(OP)/IBZ.o : $(SP)/IBZ.h $(SP)/IBZ.cpp
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/IBZ.cpp

# contains some constants and useful numerical routines
$(OP)/routines.o : $(SP)/routines.cpp $(SP)/routines.h 
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/routines.cpp

# numerical routines from NumRec
$(OP)/nrutil.o : $(SP)/nrutil.h $(SP)/nrutil.c
	$(Cpp) $(FLAGS) -c -o $@ $(SP)/nrutil.c

# clean all object and exec files
clean :
	rm -f $(RP)/$(main) $(OP)/*.o
