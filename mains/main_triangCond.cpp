#include "../source/Result.h"
#include "../source/GRID.h"
#include "../source/routines.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

int main(int argc, char *argv[])
{
  if (argc<2) exit(0);
  double T = atof(argv[1]);
  printf("T=%.3f\n",T);
  int Nkx = 1000;
  int Nnu = 1000;
  
  GRID grid("params");
  Result result(&grid);
  printf("ready\n");

  char rhoTFN[300];
  sprintf(rhoTFN,"rho.T%.3f",T);
  FILE* rhoTFile = fopen(rhoTFN,"w");
  fclose(rhoTFile);
  printf("starting\n");
  for (double U=0.6; U<4.0; U+=0.025)
  { printf("---Working: U: %.3f T:%.3f\n",U,T);
    char FN[300];
    sprintf( FN, "CHMs/CHM.U%.3f.T%.3f", U, T);
    if (not result.ReadFromFile(FN)) continue;
    double rho = result.TriangularConductivity(T, Nkx, (int)( (double)Nkx/(2.0*sqrt(3.0)) ), Nnu);
      
    FILE* rhoTFile = fopen(rhoTFN,"a");
    fprintf(rhoTFile,"%.15le %.15le %.15le\n",U,T,1.0/rho);
    fclose(rhoTFile);
  }
  return 0;
}
/*
int main()
{
  int Nkx = 1000;
  int Nnu = 800;
  
  GRID grid("params");
  Result result(&grid);

  FILE* rhoUTFile = fopen("rhoUT","w");
  fclose(rhoUTFile);
  for (double T=0.30; T>0.0099; T-=0.01)
  { for (double U=0.6; U<4.0; U+=0.025)
    { printf("---Working: U: %.3f T:%.3f\n",U,T);
      char FN[50];
      sprintf( FN, "CHM.U%.3f.T%.3f", U, T);
      result.ReadFromFile(FN);
      double rho = result.TriangularConductivity(T, Nkx, (int)( (double)Nkx/(2.0*sqrt(3.0)) ), Nnu);
      
      FILE* rhoUTFile = fopen("rhoUT","a");
      fprintf(rhoUTFile,"%.15le %.15le %.15le\n",U,T,1.0/rho);
      fclose(rhoUTFile);
    }
    FILE* rhoUTFile = fopen("rhoUT","a");
    fprintf(rhoUTFile,"\n");
  }
  fclose(rhoUTFile);
  return 0;
}
*/
