#include <cstdio>
#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"
#include <omp.h>


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


double Lorentzian(double c, double x0, double gamma, double x)
{
  return c*gamma/(sqr(x-x0)+sqr(gamma));
}

int main()
{

  GRID grid("params");
  Result result(&grid);

  grid.assign_omega(result.omega);

  Input input("params");
  double t = 0.5; // always set a default value in case parameter is not found in the input file!!!
  input.ReadParam(t,"CHM::t");

  double n=0.5;
  
  InitDelta( DOStypes::SemiCircle,	// type of DOS 
             grid.get_N(), 		// total number of omega grid points 
             0.5, 0.0, 0.01, t, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping amplitude
             result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored

  InitDOS( DOStypes::SemiCircle, 			// type of DOS 
           t, 					// hopping amplitude
           grid.get_N(), result.omega, result.NIDOS);	// total number of omega grid points, omega-grid array, the array in which NIDOS will be stored
  

  CHM chm("params");

  for(double T = 0.01; T<0.201; T+=0.01)
  for(double U = 3.5; U>2.55; U-=0.3)
  {    
    chm.SetParams(U,T,t);
    result.n=n;
    
    
    chm.Run(&result);

    char FN[50];
    sprintf( FN, "CHM.U%.3f.T%.3f", U, T );
    result.PrintResult(FN);

    char GiwFN[300];
    sprintf(GiwFN,"Giw.U%.3f.T%.3f",U,T);
    result.PrintOnImagAxis(result.G, 1000, T, GiwFN);
     

    //---------------//
    int Nw = 400;
    double w_max = 4.0;
    double* w = new double[Nw];
    for(int j=0; j<Nw; j++) w[j]=(double)j*(w_max/((double)Nw-1.0));
    double* condw = new double[Nw];

    
    //char integrandFN [300];
    //sprintf(integrandFN,"integrand.U%.3f.T%.3f.w%.3f",U,T,0.0);
    condw[0] = result.Conductivity( T , 0.5*U , 800.0, 800.0); //integrandFN );

    char DCcondUFN[300];
    sprintf(DCcondUFN,"DCcond.allomega.U%.3f",U);
    FILE* DCcondUFile = fopen(DCcondUFN,"a");
    fprintf( DCcondUFile, "%.15le %.15le %.15le\n",T,condw[0],1.0/condw[0]);
    fclose(DCcondUFile);

    printf("condw: j=0 w=0, DONE!\n");
    #pragma omp parallel for
    for(int j=1; j<Nw; j++)     
    {  //sprintf(integrandFN,"integrand.U%.3f.T%.3f.w%.3f",U,T,w[j]);
       condw[j] = result.Conductivity( w[j], T , 0.5*U , 800.0, 800.0); //integrandFN);
       printf("condw: j=%d w=%f, DONE!\n", j, w[j]);
    }

    char condwFN[300];
    sprintf(condwFN,"condw.allomega.U%.3f.T%.3f",U,T);
    PrintFunc(condwFN,Nw,condw,w);  
     

//    for(int j=0; j<Nw; j++)     
//      condw[j] = Lorentzian(0.383094,0,0.054946,w[j]);   

    char condiwFN[300];
    sprintf(condiwFN,"condiw.allomega.U%.3f.T%.3f",U,T);
    Result::PrintCondOnImagAxis(Nw, condw, w, 1000, T, condiwFN);

    

    delete [] condw;
    delete [] w;
  }

  
  return 0;
}
