#include "../source/TMT.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"
#include <cstdio>
#include <cstdlib>



// ifdef compiler directives are used to tell the compiler to ignore a block
// of code unles certain flags are used in compilation command (in this case
// -D_MPI and -D_OMP). D_OMP is used to compile the code for threaded execution
// and if both flags are used, hybrid (MPI+OpenMP) execution on multiple machines

// here, mpi.h and omp.h are not included at all if the code is compiled for serial
// execution. Otherwise, compilation would fail if those libraries are not installed
// on system.

#ifdef _MPI
#include "mpi.h"
#endif

#ifdef _OMP
#include <omp.h>
#endif 

int main(int argc, char* argv[])
{

#ifdef _MPI
          // MPI initialization. Every MPI process has an ordinal number (myrank). 
	  // 0 is the master process, and all the others slave processes
	  int myrank; 
	  MPI_Init( &argc, &argv ); 
	  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
#endif  


          // initialize grid and tmt from input file. This is done in both master and slave processes.
	  GRID grid("params");

	  TMT tmt("params");
          // NOTE: since TMT is initialized with same parameters in all process there is no actual need
	  // for master process to send the values for U and T along with Delta to slave processes.
	  // However if U and T are changed in an automatized manner (using for-loops for example) to
	  // obtain results for different sets of parameters in a single run, then this is needed since
	  // TMT in slave processes does not experience the change in parameters. Still, things can
	  // be implemented a bit differently to avoid sending U and T. For example, put the entire
	  // "if else" below in the for loop instead of just the "else" part. This way the exit signal
	  // would have to be sent after each calculation and Result would have to be initialized in both
          // master and slave processes. Maybe I'll change this in future and exclude the sending of 
          // U and T, but this way seems "cleaner" to me for now (check TMT::SolveSIAM and TMT::Slave) 


#ifdef _MPI
	  if (myrank != 0)
	  {
	    // MPI SLAVE

            // in the case of slave processes, grid has to be set maunally to the tmt 
	    // object because no result object will be passed to its Run method (neither it will be called)
	    tmt.SetGrid(&grid);
               
            // slave method sets tmt in a perpetual loop of waiting for input from master process, performing
	    // SIAM calculation on recieved data, and returning results to the master process until Exit Signal
	    // is received. If an error occurs during solving SIAM, an error signal is sent to master process.
	    tmt.Slave(myrank);
	  }
	  else 
	  {
  
	    // MPI MASTER  
#endif 
            Result result(&grid);
	    grid.assign_omega(result.omega);
          
            // create an initial guess for Delta
	    InitDelta(DOStypes::SemiCircle, grid.get_N(), 0.25, 0.0, 0.01, 0.25, result.omega, result.Delta);
            // NIDOS is not needed for input since tmt only works with bethe-specific self-consistency
    
            // opposite to what is the case with CHM, tmt takes mu as an input parameter, rather than n
            // To use tmt away from half-filling, set mu to values other than U/2
	    result.mu = tmt.get_U()/2.0;

            // this is just the initial guess for SIAM inside TMT
	    result.mu0 = 0;
             
            // run the calculation and receive the error code. 
	    bool failed = tmt.Run(&result);

	    if (failed) 
	        printf ("==== ERROR ==== Solution is INVALID\n");
        
            // prepare the file name and print out results
            char FN[50];
            sprintf(FN, "TMT.U%.3f.T%.3f.W%.3f%s", tmt.get_U(), tmt.get_T(), tmt.get_W(), 
                                                   (failed) ? ".FAILED" : "" );
	    result.PrintResult(FN);


#ifdef _MPI
            // when finished, send exit signal to slave processes to end them.
	    tmt.SendExitSignal();   
	  } //end MPI MASTER

          // all processes must reach MPI finalize.    
	  MPI_Finalize();
#endif

	  return 0;
}
