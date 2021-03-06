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

  result.n = 0.5;
  result.n0 = 0.5;
  result.mu = 0.0;
  result.mu0 = 0.0;

  for (double T=0.01; T<0.2; T+=0.01)
  {
    for (double U=0; U<4.0; U+=0.2)
    {
      InitDOS( DOStypes::Insulator, 			// type of DOS 
               t, 					// hopping amplitude
               grid.get_N(), result.omega, result.NIDOS, U);	// total number of omega grid points, omega-grid array, the array in which NIDOS will be stored

      //get_G_from_DOS(DOStypes::Insulator, t, grid.get_N(), result.omega, result.G, 0.01, U);
      //get_Sigma_from_G(t, grid.get_N(), result.omega, result.G, result.Sigma);  
      //for(int i=0; i<grid.get_N(); i++) result.Sigma = ii*0.01;

      printf("T: %.3f U: %.3f \n",U,T);

      if (T==0.01)
      {  char FN[50];
         sprintf( FN, "rigidBand.U%.3f", U );
         result.PrintResult(FN);
      }

      char intgFN[100];
      sprintf(intgFN,"integrand.U%.3f.T%.3f",U,T);

      FILE* rhocFile = fopen("rho.UT","a");
      fprintf(rhocFile, "%.15le %.15le %.15le\n", U, T, 1.0/result.NIConductivity( T , 0.0 , 600, 600, intgFN ) );
      fclose(rhocFile); 
    }
    FILE* rhocFile = fopen("rho.UT","a");
    fprintf(rhocFile, "\n");
    fclose(rhocFile);
  }

  
  return 0;
}
