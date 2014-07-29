#include <cstdio>
#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"

int main()
{
  GRID grid("params");
  Result result(&grid);
  grid.assign_omega(result.omega);
  Input input("params");
  double t = 0.5; // always set a default value in case parameter is not found in the input file!!!
  input.ReadParam(t,"CHM::t");
  
  InitDelta( DOStypes::SemiCircle,	// type of DOS 
             grid.get_N(), 		// total number of omega grid points 
             0.5, 0.0, 0.01, t, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping amplitude
             result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored

  InitDOS( DOStypes::SemiCircle, 			// type of DOS 
           t, 					// hopping amplitude
           grid.get_N(), result.omega, result.NIDOS);	// total number of omega grid points, omega-grid array, the array in which NIDOS will be stored
 

  CHM chm("params");
  
//  double U=3.0;

  for (double mu=3.0; mu<4.0; mu+=0.1)
  {
    InitDelta( DOStypes::SemiCircle,	// type of DOS 
               grid.get_N(), 		// total number of omega grid points 
               0.5, 0.0, 0.01, t, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping amplitude
               result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored

    for (double T=0.21; T<1.0; T+=0.05)
    for (double U=4.0; U<4.1; U+=0.5)
    {
      result.mu = mu;
      chm.SetParams(U,T,t);
      chm.Run(&result);
  
      char FN[50];
      sprintf( FN, "CHM.mu%.3f.U%.3f.T%.3f", result.mu, chm.get_U(), chm.get_T() );
      result.PrintResult(FN);
    }
  }

  //FILE* rhocFile = fopen("rhoc","a");
  //fprintf(rhocFile, "%.15le %.15le\n", chm.get_T(), result.Conductivity( chm.get_T() , 0.5*chm.get_U() , 400.0, 400.0, 10000.0 ) );
  //fclose(rhocFile);
  //PrintFunc("Sigma", grid.get_N(), result.Sigma, result.omega);
  
  return 0;
}
