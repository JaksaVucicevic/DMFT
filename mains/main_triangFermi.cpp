#include "../source/Result.h"
#include "../source/GRID.h"
#include "../source/routines.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

double epsilon(double kx, double ky)
{
  return -2.0 * 0.25 * ( cos(kx) + 2.0 * cos(kx/2.0) * cos(sqrt(3.0)*ky/2.0) ); 
}

int main()
{ 
  GRID grid("params");
  Result result(&grid);

  for (double T=0.02; T>0.005; T-=0.06)
  { for (double U=3.9; U>3.5; U-=0.1)
    { printf("---Working: U: %.3f T:%.3f\n",U,T);
      char FN[50];
      sprintf( FN, "CHM.U%.3f.T%.3f", U, T);
      if (not result.ReadFromFile(FN)) continue;

      char AFN[50];
      sprintf( AFN, "Fermi.U%.3f.T%.3f", U, T);
      result.PrintFermiSurface(AFN, &epsilon);
    }
  }
  return 0;
}
