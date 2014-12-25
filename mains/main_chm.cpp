#include <cstdio>
#include "../source/LambdaCalculator.h"
#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"

int main()
{
  // initialize grid with params read from the input file "params"
  GRID grid("params");

  // initialize result object for storing all freq. dependent functions on the above grid
  Result result(&grid);

  // fill in the omega grid of the object
  grid.assign_omega(result.omega);

  // read the hopping amplitude that will be used in CHM from the input file, to use it in initializeing Delta and NIDOS
  Input input("params");
  double t = 0.5; // always set a default value in case parameter is not found in the input file!!!
  input.ReadParam(t,"CHM::t");
  

  // initialize hybridization function that will be used for the initial guess 
  InitDelta( DOStypes::SemiCircle,	// type of DOS 
             grid.get_N(), 		// total number of omega grid points 
             t, 0.0, 0.01, t, 	// hybridization-V, shift, broadening eta, and hopping amplitude (used in DOS)
             result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored

  // initialize the non-interacting DOS (NIDOS) that will be used for calculating G 
  // (...unless bethe-specific self-consistency is used. in that case you don't need this function call)
  InitDOS( DOStypes::SemiCircle, 			// type of DOS 
           t, 					// hopping amplitude
           grid.get_N(), result.omega, result.NIDOS);	// total number of omega grid points, omega-grid array, the array in which NIDOS will be stored

 
  // create a CHM object with params set in the input file "params". Those params not found in the input file will be set to their default values
  CHM chm("params");
  chm.LC->SetDoOutput(false);

  //iterate over parameters
  
  for (double T=0.1; T<0.3; T+=0.02)
  { 
    //reset the file in which resistivity and lambda will be printed
    char rhoTFN[300];
    sprintf(rhoTFN,"rho.T%.3f",T);
    FILE* rhoTFile = fopen(rhoTFN,"w");
    fclose(rhoTFile);

    for (double U=2.0; U<4.0; U+=0.2)   
    //  for (double n=0.5; n<0.9; n+=0.02)
    {  double n = 0.5;
       result.n = n;  
       //initialize Delta again if you do not want to start from the previous solution
/*     InitDelta( DOStypes::SemiCircle,	// type of DOS 
                  grid.get_N(), 		// total number of omega grid points 
                  0.5, 0.0, 0.01, t, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping amplitude
                  result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored
*/

       // set parameters of the calculation
       chm.SetParams(U,T,t);
       // perform the DMFT CHM calculation on the result object in which initial Delta, NIDOS and n are now found and will be used as input
       chm.Run(&result);

       //get lambda
       double best_lambda = chm.LC->best_lambda;

       // now print the result to a file named accordingly to the parameters of the calculation
       char FN[300];
       sprintf( FN, "CHM.n%.3f.U%.3f.T%.3f", n, U, T );
       result.PrintResult(FN);

       // print conductivity and lambda to a file for given T
       FILE* rhoTFile = fopen(rhoTFN,"a");
       fprintf(rhoTFile, "%.15le %.15le %.15le\n", U, 1.0/result.Conductivity(T, 0.5*U, 400, 400, NULL, false), best_lambda );
       fclose(rhoTFile);
    }
  }
  
  return 0;
}
