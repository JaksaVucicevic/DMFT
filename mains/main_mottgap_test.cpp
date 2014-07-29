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

  
  for (double U=2.6; U<3.5; U+=0.025)
  {
    char rhoFN[300];
    sprintf(rhoFN,"test_rho_vs_T.U%.3f",U);
    FILE* rhoFile = fopen(rhoFN,"w");
 
    double ltco = 1.0/50.0;
    if (U > 2.949) ltco=1.0/35.0;
    if (U > 3.124) ltco=1.0/32.0;
    if (U > 3.249) ltco=1.0/30.0;
    if (U > 3.399) ltco=1.0/27.0;

    bool found = false;
    for (double T=0.005; T<0.25; T+=0.005)
    { if (T<ltco) continue;

      if (not found)        
      {  char FN[300];
         sprintf( FN, "CHM.U%.3f.T%.3f", U, T );
         result.ReadFromFile(FN);

         char testFN[300];
         sprintf( testFN, "test_CHM.U%.3f.T%.3f", U, T );       
         result.PrintResult(testFN);
         
      }
      found = true;

      fprintf(rhoFile, "%.15le %.15le %.15le\n", T,
                   result.Conductivity( T , U*0.5 , 400, 400, NULL, false ),
                   result.Conductivity( T , U*0.5 , 400, 400, NULL, true )
              );
    }
    fclose(rhoFile);
  }

  
  return 0;
}
