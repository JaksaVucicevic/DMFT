#include <cstdio>
#include <cstdlib>
#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"

int main(int argc, char* argv[])
{
  if (argc!=2) exit(0);
  double T = atof(argv[1]);
  printf(" T=%.3f\n",T);
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
  //for (double T=0.30; T>0.0099; T-=0.01)
  { 
    result.mu = 0.25;
  
    InitDelta( DOStypes::Uniform,	// type of DOS 
               grid.get_N(), 		// total number of omega grid points 
               0.5, 0.375, 0.05, 1.125, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping
               result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored
    printf(" initialized Delta...\n");
    //result.ReadFromFile("CHM.U2.500.T0.100");
    //result.ReadFromFile("CHM.U0.500.T0.240");
    result.n=n;
    result.mu0 = 0.0;
    chm.SetSIAMeta(0);
    chm.SetMixerOptions(2, (int []){1,0});

    bool alreadyBroken = false;
    result.mu=1.403676e+00;
    for (double U=3.0; U<5.0; U+=0.500)
    {
      Result rcpy(result);
      chm.SetParams(U,T,0);
      bool not_conv = chm.Run(&result);
      /*if ((not_conv)and(not alreadyBroken))
      { alreadyBroken = true;
        result.CopyFrom(rcpy);
        chm.SetSIAMeta(1e-2, true);
        chm.SetMixerOptions(4, (int []){1,2,2,1});
        not_conv = chm.Run(&result);
      }*/
      char FN[50];
      sprintf( FN, "CHM.U%.3f.T%.3f%s", chm.get_U(), chm.get_T(), (not_conv) ? ".FAILED" : "" );
      result.PrintResult(FN);
      //if (not_conv) break;
    }
  } 
  return 0;
}
