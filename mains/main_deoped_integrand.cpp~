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

//  double n = 0.5; // always set a default value in case parameter is not found in the input file!!!
//  input.ReadParam(n,"main::n");

  double U=1.0;
  for(double T=0.20; T>=0.199; T-=0.02)
  {
  for(double n=0.960; n<0.999; n+=0.002)
  {
    result.n = n;
  
    CHM chm("params");
    result.n = n;
    chm.SetParams(U,T,t);

    //chm.Run(&result);
  
    // now print the result to a file named accordingly to the parameters of the calculation
    char FN[50];
    sprintf( FN, "CHM.n%.3f.U%.3f.T%.3f", n, chm.get_U(), chm.get_T() );
    //result.PrintResult(FN);
    result.ReadFromFile(FN);
    sprintf( FN, "integrand.n%.3f.U%.3f.T%.3f", n, chm.get_U(), chm.get_T() );
    result.Conductivity( chm.get_T(), result.mu , 400.0, 400.0, 10000.0, FN );
    //FILE* rhocFile = fopen("rho.nT.large_n","a");
    //fprintf(rhocFile, "%.15le %.15le %.15le %.15le\n",  result.mu, n, chm.get_T(), result.Conductivity( chm.get_T(), result.mu , 400.0, 400.0, 10000.0 ) );
    //fclose(rhocFile);
    
  }
  //FILE* rhocFile = fopen("rho.nT.large_n","a");
  //fprintf(rhocFile, "\n");
  //fclose(rhocFile);
    
  }
  return 0;
}
