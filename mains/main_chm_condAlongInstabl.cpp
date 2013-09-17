#include <cstdio>
#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"
#include <omp.h>
/*
double Us [] = {
2.4236347826,
2.3896000000,
2.3573565217,
2.3322782609,
2.3143652174,
2.2946608696,
2.2803304348,
2.2660000000,
2.2570434783,
2.2552521739,
2.2552521739,
2.2624173913,
2.2749565217,
2.2964521739,
2.3251130435,
2.3555652174,
2.3860173913,
2.4200521739,
2.4558782609,
2.4952869565,
2.5418608696};

double Ts [] = {
0.0519810268,
0.0566015625,
0.0619151786,
0.0663046875,
0.0711562500,
0.0764698661,
0.0817834821,
0.0884832589,
0.0956450893,
0.1018828125,
0.1090446429,
0.1171305804,
0.1261406250,
0.1351506696,
0.1448537946,
0.1547879464,
0.1640290179,
0.1730390625,
0.1827421875,
0.1910591518,
0.2005312500};
*/

double Us [] = {
2.31822,
2.25683,
2.2657,
2.31034,
2.37233,
2.44513};

double Ts [] = {
0.07,
0.096,
0.12,
0.14,
0.16,
0.18};

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

  for(int i = sizeof(Us)/sizeof(double) - 1; i>=0; i--)
  //for(double dU = -0.6; dU<0.61; dU+=0.3)
  { double dU = 0.0;
    double n = 0.5; // always set a default value in case parameter is not found in the input file!!!
    input.ReadParam(n,"main::n");
    result.n = n;

    chm.SetParams(Us[i]+dU,Ts[i],t);
    chm.Run(&result);
  
    char FN[50];
    sprintf( FN, "CHM.n%.3f.T%.3f.dU%.3f", n, chm.get_T(), dU );
    result.PrintResult(FN);

    //---------------//
/*    int Nw = 51;
    double dw = 0.1;
    double* condw = new double[Nw];
    condw[0] = result.Conductivity( chm.get_T() , 0.5*chm.get_U() , 400.0, 400.0, 10000.0 );
    #pragma omp parallel for
    for(int j=1; j<Nw; j++)     
    {  condw[j] = result.Conductivity( j*dw, chm.get_T() , 0.5*chm.get_U() , 400.0, 400.0, 10000.0 );
       printf("condw: j=%d w=%f, DONE!\n", j, j*dw);
    }

    for(int j=0; j<Nw; j++)     
    {  FILE* rhocFile = fopen("rhow","a");
       fprintf(rhocFile, "%d %.15le %.15le %.15le\n", i, Ts[i], j*dw, condw[j] );
       fclose(rhocFile);     
    }
    FILE* rhocFile = fopen("rhow","a");
    fprintf(rhocFile, "\n");
    fclose(rhocFile);    
 
    delete [] condw;
*/
    //---------------//
/*    char SigmasFN[100];
    sprintf(SigmasFN,"Sigmas.dU%.3f",dU);
    char DOSesFN[100];
    sprintf(DOSesFN,"DOSes.dU%.3f",dU);

    for(int j=0; j<grid.get_N(); j++)
    {  FILE* SigmasFile = fopen(SigmasFN,"a");
       fprintf(SigmasFile, "%d %.15le %.15le %.15le %.15le\n", i, Ts[i], result.omega[j], real(result.SOCSigma[j]), imag(result.Sigma[j]) );
       fclose(SigmasFile);
       FILE* DOSesFile = fopen(DOSesFN,"a");
       fprintf(DOSesFile, "%d %.15le %.15le %.15le\n", i, Ts[i], result.omega[j], -(1.0/3.1415)*imag(result.G[j]) );
       fclose(DOSesFile);          
    }
    FILE* SigmasFile = fopen(SigmasFN,"a");
    fprintf(SigmasFile, "\n");
    fclose(SigmasFile);       
    FILE* DOSesFile = fopen(DOSesFN,"a");
    fprintf(DOSesFile, "\n");
    fclose(DOSesFile);     
*/
    //---------------//

  }

  
  return 0;
}
