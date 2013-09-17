#include "/home/jaksa/DMFT/source/Result.h"
#include "/home/jaksa/DMFT/source/GRID.h"
#include "/home/jaksa/DMFT/source/routines.h"
#include <cstdio>

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
 0.0084365782
};


int main()
{
  // initialize grid. params are read from input file
  GRID grid("params");

  // initialize a result object on grid


  FILE* A0UTFile = fopen("ImGiw1","w");
//  for(double U=2.4; U<3.501; U+=0.05)
//  for(double T=0.002; T<0.08; T+=0.002)
  for(int i=0; i<sizeof(Us)/sizeof(double); i++)
//  for(int i=sizeof(Us)/sizeof(double)-1; i>=0; i--)
  {
    double U = Us[i];
    double T = Ts[i];  
    
    printf("Working: U: %.3f T: %.5f\n",U,T);
    char FN[200];
    sprintf(FN,"CHM.n0.500.U%.3f.T%.5f.Met",U,T);
    Result result(&grid);
    result.ReadFromFile(FN);
    
    sprintf(FN,"IAGf.U%.3f.T%.3f",U,T);
    result.PrintOnImagAxis(result.G, 1024, T, FN);       

    double ImGw0 = imag(result.G[grid.get_N()/2]);

    double iw1, ReGiw1, ImGiw1;

    FILE* f = fopen(FN,"r"); 
    fscanf(f,"%le %le %le", &iw1, &ReGiw1, &ImGiw1);
    fclose(f);

    fprintf(A0UTFile, "%le %le %le %le\n",U,T, ImGw0, ImGiw1);

    
    //sprintf(FN,"ImagAxis2/IASigf.U%.3f.T%.3f",U,T);
    //result.PrintOnImagAxis(result.Sigma, 4096, T, FN);     

  }
  fclose(A0UTFile);
  return 0;
}
