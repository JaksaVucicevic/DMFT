#include <cstdio>
#include <complex>
#include "../source/routines.h"
#include "../source/IBZ.h"
#include "../source/GRID.h"


int main()
{

  int N=1000;
  GRID grid(N, 1.5, false);
  
  //GRID grid("params");
  //int N=grid.get_N();

  double* omega = new double[N];
  double* dos = new double[N];

  grid.assign_omega(omega);
  printf("-----grid ready\n");

  int Nx = 4000;
  int Ny = 4000;
  IBZ ibz(IBZtypes::TriangularLattice, Nx, Ny);
  printf("-----ibz ready\n");
  for(int n=0; n<N; n++)
  { for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++)
      ibz.summand[i][j] = 1.0/(omega[n]-ibz.epsilon[i][j]+ii*1e-3);
    printf("-----summand ready, n=%d \n",n);
    dos[n] = -(1.0/pi)*imag(ibz.sum());    
  }
  PrintFunc("triangular_dos",N,dos,omega);
  printf("Spectral weight: %.6f\n",TrapezIntegral(N, dos, omega));
  delete [] dos;  
  delete [] omega;
  return 0;
}


/*
int main()
{
  int N=500;
  double* omega = new double[N];
  double* dos = new double[N];
  GRID grid(N, 2.0, false);
  grid.assign_omega(omega);
  printf("-----grid ready\n");
  int Nx = 4000;
  int Ny = 4000;
  IBZ ibz(IBZtypes::TriangularLattice, Nx, Ny);
  printf("-----ibz ready\n");
  for(int n=0; n<N; n++)
  { for(int i=0; i<Nx; i++)
    for(int j=0; j<Ny; j++)
      ibz.summand[i][j] = 1.0/(omega[n]-ibz.epsilon[i][j]+ii*1e-3);
    printf("-----summand ready, n=%d \n",n);
    dos[n] = -(1.0/pi)*imag(ibz.sum());    
    if (n%20==0)
    {
      char FN[300];
      sprintf(FN,"Ak.w%.3f",omega[n]);
      printf("-----print ready ready, n=%d \n",n);
      ibz.PrintToFile(FN);
    }
  }
  PrintFunc("dos",N,dos,omega);
  printf("Spectral weight: %.6f\n",TrapezIntegral(N, dos, omega));
  delete [] dos;  
  delete [] omega;
  return 0;
}
*/
