#include "../source/TMT.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"
#include "../source/LambdaCalculator.h"

#include <cstdio>
#include <cstdlib>

#ifdef _MPI
#include "mpi.h"
#endif

#ifdef _OMP
#include <omp.h>
#endif 

int main(int argc, char* argv[])
{

#ifdef _MPI
	  int myrank; 
	  MPI_Init( &argc, &argv ); 
	  MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
#endif  

	  GRID grid("params");

	  TMT tmt("params");

#ifdef _MPI
	  if (myrank != 0)
	  {
	    // MPI SLAVE

	    tmt.SetGrid(&grid);
               
	    tmt.Slave(myrank);
	  }
	  else 
	  {
  
	    // MPI MASTER  
#endif 
            Result result(&grid);
	    grid.assign_omega(result.omega);
            
            InitDOS(DOStypes::SemiCircle, 0.5, grid.get_N(), result.omega, result.NIDOS);

            for(double T=0.05; T<0.06; T += 10000.0)
            for(double U=0.25; U<=7.1; U+=0.5)
            {
	      InitDelta(DOStypes::SemiCircle, grid.get_N(), 0.5, 0.0, 0.01, 0.5, result.omega, result.Delta);
              
              for(double W=0.25; W<=7.1; W+=0.5)
              {   
                //InitDelta(DOStypes::SemiCircle, grid.get_N(), 0.5, 0.0, 0.01, 0.5, result.omega, result.Delta);
                
                result.mu = U/2.0;
            
                result.mu0 = 0;

                tmt.SetWDN(W, Distributions::Uniform, (W>0.0) ? 20 + (int) ( ((U>2.8) ? 2.0 : 1.0) * W / 0.1 ) : 1);
                tmt.SetParams(U, T, 0.5);

                char OFN[300];
                sprintf(OFN,".T%.3f.U%.3f.W%.3f",T,U,W);             
                tmt.LC->SetOutputFileName(OFN);
                tmt.LC->SetNMT(grid.get_N(), 1000, T);
                //printf("grid.get_N(): %d, half: %d\n",grid.get_N(),grid.get_N()/2) ;
                //tmt.LC->SetOffset(grid.get_N()/2);

                bool failed = tmt.Run(&result);

                FILE* f = fopen("lambdaUW.T0.050","a");
                fprintf(f, "%.15le %.15le %.15le %.15le\n", U,W, tmt.LC->lambdas[0], tmt.LC->continued_lambdas[0]);
                fclose(f);
                 
                if (failed) 
                  printf ("==== ERROR ==== Solution is INVALID\n");
        
                char FN[50];
                sprintf(FN, "TMT.T%.3f.U%.3f.W%.3f%s", tmt.get_T(), tmt.get_U(), tmt.get_W(), 
                                                       (failed) ? ".FAILED" : "" );
	        result.PrintResult(FN);
             }
           }


#ifdef _MPI
	    tmt.SendExitSignal();   
	  } //end MPI MASTER

	  MPI_Finalize();
#endif

	  return 0;
}
