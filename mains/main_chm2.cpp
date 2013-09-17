#include <cstdio>
#include "../source/SIAM2.h"
#include "../source/CHM2.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"

int main()
{
  // initialize grid with params read from the input file "params"
  GRID grid("params");

  // initialize result object for discretizing relevant functions on the above grid
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
             0.5, 0.0, 0.01, t, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping amplitude
             result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored

  // initiali the non-interacting DOS (NIDOS) that will be used for calculating G 
  // (...unless bethe-specific self-consistency is used. in that case you don't need this function call)
  InitDOS( DOStypes::SemiCircle, 			// type of DOS 
           t, 					// hopping amplitude
           grid.get_N(), result.omega, result.NIDOS);	// total number of omega grid points, omega-grid array, the array in which NIDOS will be stored

  // read the occupation number and put it in result to be used as input for CHM
  double n = 0.5; // always set a default value in case parameter is not found in the input file!!!
  input.ReadParam(n,"main::n");
  result.n = n;

  // create a CHM object with params set in the input file "params". Those params not found in the input file will be set to their default values
  CHM2 chm("params");

//  double Us [] = {2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0};
//  double mus [] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6};

  double Us [] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.4, 3.5, 3.7};
  double mus [] = {2.0, 2.01, 2.02, 2.03, 2.04};
  //double mus [] = {1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0};

  result.n = 0.5;
  result.n0 = 0.5;
  
  for(int i=1; i<5; i++)
  {
    double U = 4.0;
    double mu = mus[i];
    double T = 0.02;
    result.mu = mu;
    result.mu0 = mu;

    chm.SetParams(U,T,t);
    if (i==0)
    {  chm.siam->HalfFilling = true;
       chm.SetForceSymmetry(true);
    }
    else
    {  chm.siam->HalfFilling = false;
       chm.SetForceSymmetry(false);
    }
    
    // perform the DMFT CHM calculation on the result object in which initial Delta, NIDOS and n are now found and will be used as input
    chm.Run(&result);
  
  // now print the result to a file named accordingly to the parameters of the calculation
    char FN[50];
    sprintf( FN, "CHM.U%.3f.T%.3f.mu%.3f.n%.3f", U, T, mu, result.n );
    result.PrintResult(FN);
  }

  //FILE* rhocFile = fopen("rhoc","a");
  //fprintf(rhocFile, "%.15le %.15le\n", chm.get_T(), result.Conductivity( chm.get_T() , 0.5*chm.get_U() , 400.0, 400.0, 10000.0 ) );
  //fclose(rhocFile);
  //PrintFunc("Sigma", grid.get_N(), result.Sigma, result.omega);
  
  return 0;
}
