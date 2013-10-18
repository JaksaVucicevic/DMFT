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

            for(double T=0.07; T>0.009; T -= 0.02)
            for(double U=1.5; U<3.6; U+=0.5)
            {
	      InitDelta(DOStypes::SemiCircle, grid.get_N(), 0.5, 0.0, 0.01, 0.5, result.omega, result.Delta);
              
              for(double W=0.0; W<4.1; W+=0.2)
              {   
                //InitDelta(DOStypes::SemiCircle, grid.get_N(), 0.5, 0.0, 0.01, 0.5, result.omega, result.Delta);
                
                result.mu = U/2.0;
            
                result.mu0 = 0;

                tmt.SetWDN(W, Distributions::Uniform, (W>0.0) ? (int) (W/0.1) : 1);
                tmt.SetParams(U, T, 0.5);

                char OFN[300];
                sprintf(OFN,".U%.3f.T%.3f.W%.3f",U,T,W);             
                tmt.LC->SetOutputFileName(OFN);
                tmt.LC->SetT(T);
                bool failed = tmt.Run(&result);

                FILE* f = fopen("lambda.UTW","a");
                fprintf(f, "%.15le %.15le %.15le %.15le\n",tmt.LC->continued_lambdas[0]);
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
