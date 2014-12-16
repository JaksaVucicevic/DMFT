#include <cstdio>
#include <cstdlib>
#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"

int main(int argc, char* argv[])
{
  /*if (argc!=2) exit(0);
  double T = atof(argv[1]);
  printf(" T=%.3f\n",T);*/
  GRID grid("params");
  Result result(&grid);
  grid.assign_omega(result.omega);

  Input input("params");
  ReadFunc("triangular_dos", grid.get_N(), result.NIDOS, result.omega);

  // read the occupation number and put it in result to be used as input for CHM
  double n = 0.5+1e-10; // always set a default value in case parameter is not found in the input file!!!
  result.n=n;
  result.mu0 = 0.0;
  printf(" about to start...\n");
  // create a CHM object with params set in the input file "params". Those params not found in the input file will be set to their default values
  CHM chm("params");
  for (double T=0.18; T>0.05; T-=0.03)
  { chm.UseFixedMuSIAMRun=false;
    result.mu = 0.25;
    result.n=n;
    //result.mu0 = 0.0;
    chm.SetSIAMeta(0, false);
    chm.SetMixerOptions(2, (int []){1,0});

    bool alreadyBroken = false;

    for (double U=1.5; U<3.0; U+=0.1)
    {
      InitDelta( DOStypes::Uniform,	// type of DOS 
                 grid.get_N(), 		// total number of omega grid points 
                 0.5, 0.375, 0.05, 1.125, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping
                 result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored

      Result rcpy(result);
      chm.SetParams(U,T,0);
      bool not_conv = chm.Run(&result);
 
      char FN[50];
      sprintf( FN, "CHM.U%.3f.T%.3f%s", chm.get_U(), chm.get_T(), (not_conv) ? ".FAILED" : "" );
      result.PrintResult(FN);
      //if (not_conv) break;
    }
  } 
  return 0;
}
