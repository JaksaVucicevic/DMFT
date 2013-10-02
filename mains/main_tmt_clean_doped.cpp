#include "../source/TMT.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"
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

            //for(double T=0.02; T>0.009; T -= 0.02)
            double T=0.02;
            double U=4.0;
            //for(double U=3.0; U<4.1; U+=0.5)
            {
	      InitDelta(DOStypes::SemiCircle, grid.get_N(), 0.5, 0.0, 0.01, 0.5, result.omega, result.Delta);
            
              //for(double W=0.6; W<4.1; W+=0.2)
              double W=0;
              for(double dmu=0.0; dmu<2.0; dmu+=1000.05)
              {   
                //InitDelta(DOStypes::SemiCircle, grid.get_N(), 0.5, 0.0, 0.01, 0.5, result.omega, result.Delta);
                
                result.mu = U/2.0+dmu;
            
                result.mu0 = 0;

                tmt.SetWDN(W, Distributions::Uniform, (W>0.0) ? (int) (W/0.05) : 1);
                tmt.SetParams(U, T, 0.5);             
                bool failed = tmt.Run(&result);

                if (failed) 
                  printf ("==== ERROR ==== Solution is INVALID\n");
        
                char FN[50];
                sprintf(FN, "TMT.T%.3f.U%.3f.dmu%.3f%s", tmt.get_T(), tmt.get_U(), dmu, 
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
