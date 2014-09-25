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
  

  InitDOS( DOStypes::SemiCircle, 			// type of DOS 
           t, 					// hopping amplitude
           grid.get_N(), result.omega, result.NIDOS);	// total number of omega grid points, omega-grid array, the array in which NIDOS will be stored

  double n = 0.5; // always set a default value in case parameter is not found in the input file!!!
  result.n = n;
  
  CHM chm("params");
  for (double U=2.6; U<3.5; U+=0.025)
  {
    char rhoFN[300];
    sprintf(rhoFN,"UHTMI_rho_vs_T.U%.3f",U);
    FILE* rhoFile = fopen(rhoFN,"w");

    InitDelta( DOStypes::SemiCircle,	// type of DOS 
               grid.get_N(), 		// total number of omega grid points 
               0.5, 0.0, 0.01, t, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping amplitude
               result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored

    for (double T=1.0; T>0.05; T-=0.05)
    {
      result.n = n;
      chm.SetParams(U,T,t);

      chm.Run(&result);
  
      char FN[300];
      sprintf( FN, "CHM.U%.3f.T%.3f", chm.get_U(), chm.get_T() );
      result.PrintResult(FN);

      fprintf(rhoFile, "%.15le %.15le %.15le\n", T,
                   result.Conductivity( T , U*0.5 , 400, 400, NULL, false ),
                   result.Conductivity( T , U*0.5 , 400, 400, NULL, true )
              );
    }
    fclose(rhoFile);
  }

  
  return 0;
}
