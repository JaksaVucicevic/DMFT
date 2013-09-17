#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"

int main()
{

  double Ts [] = 
  {
	7.000000000000001e-02, 
	8.000000000000000e-02, 
	9.000000000000000e-02, 
	9.999999999999999e-02, 
	1.100000000000000e-01, 
	1.200000000000000e-01, 
	1.300000000000000e-01, 
	1.400000000000000e-01, 
	1.500000000000000e-01, 
	1.600000000000000e-01, 
	1.700000000000000e-01, 
	1.800000000000000e-01, 
	1.900000000000000e-01 
  };

  double Us [] = 
  {
	2.359999999999992,
	2.329999999999993,
	2.299999999999994,
	2.279999999999994,
	2.259999999999994,
	2.239999999999995,
	2.219999999999995,
	2.199999999999996,
	2.179999999999996,
	2.169999999999996,
	2.149999999999997,
	2.129999999999997,
	2.119999999999997
  };
  int NT = sizeof(Ts)/sizeof(double);


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

  double n = 0.5; // always set a default value in case parameter is not found in the input file!!!
  input.ReadParam(n,"main::n");
  result.n = n; 

  CHM chm("params");


  for (int i=0; i<NT; i++)
  {
    chm.SetParams(Us[i],Ts[i],t);
    chm.Run(&result);
  
    char FN[50];
    sprintf( FN, "CHM.U%.3f.T%.3f", chm.get_U(), chm.get_T() );
    result.PrintResult(FN);

    FILE* rhocFile = fopen("rhoc","a");
    double rhoc = result.Conductivity( chm.get_T() , 0.5*chm.get_U() , 400.0, 400.0, 10000.0 ) ;
    fprintf(rhocFile, "%.15le %.15le %.15le\n", chm.get_T(),  rhoc, 1/rhoc);
    fclose(rhocFile);

    for (double dU=-0.3; dU<0.31; dU+=0.04)
    {
      chm.SetParams(Us[i]+dU,Ts[i],t);
      chm.Run(&result);
  
      char FN[50];
      sprintf( FN, "CHM.U%.3f.T%.3f", chm.get_U(), chm.get_T() );
      result.PrintResult(FN);

      char rhoFN[50];
      sprintf( rhoFN, "rho.dU%.3f", dU );
      
      FILE* rhoFile = fopen(rhoFN,"a");
      double rho = result.Conductivity( chm.get_T() , 0.5*chm.get_U() , 400.0, 400.0, 10000.0 ) ;
      fprintf(rhoFile, "%.15le %.15le %.15le %.15le %.15le\n", chm.get_T(),  rho, 1.0/rho, rho/rhoc, rhoc/rho);
      fclose(rhoFile);
    }
  }

  
  return 0;
}
