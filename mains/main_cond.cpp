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
  

  for (double U=3.0; U<4.1; U+=0.5)
  { char condnTFN[50];
    sprintf( condnTFN, "cond.U%.3f.nT", U);
    FILE* condnTFile = fopen(condnTFN,"w");

    for (double n=0.525; n<0.75; n+=0.025)
    {
      char condnFN[50];
      sprintf( condnFN, "cond.U%.3f.n%.3f", U, n );
      FILE* condnFile = fopen(condnFN,"w");

      for (double T=0.01; T<1.0; T+=0.05)
      {
        char FN[50];
        sprintf( FN, "CHM.n%.3f.U%.3f.T%.3f", n, U, T );
        result.ReadFromFile(FN);
        printf( "%s ----- mu = %.3f mu0 = %.3f", FN, result.mu, result.mu0);

        char integrandFN[50];
        sprintf( integrandFN, "integrands/integrand.n%.3f.U%.3f.T%.3f", n, U, T );
        double sigma = result.Conductivity(T, result.mu, 400, 400, integrandFN);

        fprintf(condnTFile,"%.15le %.15le %.15le %.15le\n", n, T, 1.0/sigma, result.mu);
        fprintf(condnFile,"%.15le %.15le %.15le\n", T, 1.0/sigma, result.mu);
      }
      fclose(condnFile);
      fprintf(condnTFile,"\n");
    }
    fclose(condnTFile);
  }
  
  //FILE* rhocFile = fopen("rhoc","a");
  //fprintf(rhocFile, "%.15le %.15le\n", chm.get_T(), result.Conductivity( chm.get_T() , 0.5*chm.get_U() , 400.0, 400.0, 10000.0 ) );
  //fclose(rhocFile);
  //PrintFunc("Sigma", grid.get_N(), result.Sigma, result.omega);
  
  return 0;
}
