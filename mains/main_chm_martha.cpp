#include <cstdio>
#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"

double Us [] = {
/*2.3876262626,
2.4166666667,
2.4318181818,
2.4469696970,
2.4621212121,
2.4734848485,
2.4823232323,
2.4936868687,
2.5025252525,
2.5113636364,
2.5214646465,
2.5315656566,
2.5454545455,
2.5606060606,
2.5606060606,
2.5909090909,
2.6300505051,
2.6679292929,
2.7070707071,
2.7550505051,
2.8005050505,
2.8472222222,*/
2.8901515152,
2.9330808081,
2.9760101010,
3.0239898990,
3.0694444444,
3.1047979798,
3.1578282828,
3.1994949495,
3.2323232323
//4.0
};
double Ts [] = {
/* 0.0522123894,
 0.0502064897,
 0.0489675516,
 0.0479056047,
 0.0465486726,
 0.0456637168,
 0.0448967552,
 0.0440707965,
 0.0432448378,
 0.0424188791,
 0.0416519174,
 0.0407669617,
 0.0398230088,
 0.0385250737,
 0.0385250737,
 0.0364601770,
 0.0343362832,
 0.0323303835,
 0.0302064897,
 0.0278466077,
 0.0256637168,
 0.0235988201,*/
 0.0216519174,
 0.0198230088,
 0.0180530973,
 0.0163421829,
 0.0146312684,
 0.0132153392,
 0.0113274336,
 0.0097345133,
 0.0084365782,
// 0.01
};

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

  double n = 0.5; // always set a default value in case parameter is not found in the input file!!!
  input.ReadParam(n,"main::n");
  result.n = n;
  
  FILE* Aw0MetFile = fopen ("Aw0_MET","w");
  fclose(Aw0MetFile);

  for(int i=0; i<sizeof(Us)/sizeof(double); i++)
//  for(int i=sizeof(Us)/sizeof(double)-1; i>=0; i--)
  {
    double U = Us[i];
    double T = Ts[i];  
   
    CHM chm("params");
    result.n = n;
    chm.SetParams(U,T,t);

    bool err = chm.Run(&result);
  
    FILE* Aw0MetFile = fopen ("Aw0_MET","a");
    fprintf(Aw0MetFile,"%.15le %.15le %.15le %.15le %.15le\n", U,T,imag(result.G[grid.get_N()/2]),-2.0*imag(result.Sigma[grid.get_N()/2]), result.get_ImGiw1(T));
    fclose(Aw0MetFile);
    
    char FN[50];
    sprintf( FN, "CHM.n%.3f.U%.3f.T%.5f.Met", n, chm.get_U(), chm.get_T() );
    result.PrintResult(FN);

    if (err) 
      InitDOS( DOStypes::SemiCircle, 			// type of DOS 
               t, 					// hopping amplitude
               grid.get_N(), result.omega, result.NIDOS);	// total number of omega grid points, omega-grid array, the array in which NIDOS will be stored
  }

  FILE* Aw0InsFile = fopen ("Aw0_INS","w");
  fclose(Aw0InsFile);


//  for(int i=0; i<sizeof(Us)/sizeof(double); i++)
  for(int i=sizeof(Us)/sizeof(double)-1; i>=0; i--)
  {
    double U = Us[i];
    double T = Ts[i];  
   
    CHM chm("params");
    result.n = n;
    chm.SetParams(U,T,t);

    bool err = chm.Run(&result);
  
    FILE* Aw0InsFile = fopen ("Aw0_INS","a");
    fprintf(Aw0InsFile,"%.15le %.15le %.15le %.15le %.15le\n", U,T,imag(result.G[grid.get_N()/2]),-2.0*imag(result.Sigma[grid.get_N()/2]), result.get_ImGiw1(T));
    fclose(Aw0InsFile);
    
    char FN[50];
    sprintf( FN, "CHM.n%.3f.U%.3f.T%.5f.Ins", n, chm.get_U(), chm.get_T() );
    result.PrintResult(FN);

    if (err) 
      result.ReadFromFile("CHM.n0.500.U4.000.T0.01000.Ins") ;

  }
  
  return 0;
}
