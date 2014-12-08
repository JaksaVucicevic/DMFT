#include <cstdio>
#include <cmath>
#include "../source/routines.h"
#include "../source/GRID.h"
/*
const int N = 230;

double f(double x)
{
  return DOS(DOStypes::SquareLattice, 0.5, x, 0.0);
}

int main()
{
  double shift = 0.0;

  double* X = new double[N];
  GRID grid(0, N, 5.0, 0.0, 0.0);
  grid.assign_omega(X);
  
  double* Y = new double[N];
  InitDOS(DOStypes::SquareLattice, 0.5, N, X, Y);
  printf("TrapezIntegral: %.6f\n", TrapezIntegral(N, Y, X) ); 


  for(int i=0; i<N; i++) X[i] += shift;
  
  printf("SmartIntegral: %.6f\n", SmartIntegral(N, Y, X, shift, 0.5, &f, 300, 1e-8) ); 

  return 0;
}
*/

/*

const int N = 230;
double x0=4.8;

double f(double x)
{
  return 1.0/(exp(x-x0)-1.0);
}

int main()
{
  double shift = 0.0;

  double* X = new double[N];
  GRID grid(0, N, 5.0, 0.0, 0.0);
  grid.assign_omega(X);
  
  double* Y = new double[N];
  for(int i=0; i<N; i++) Y[i] = f(X[i]);
  PrintFunc("raw",N,Y,X);
  printf("TrapezIntegral: %.6f\n", TrapezIntegral(N, Y, X) ); 
 
  printf("SmartIntegral: %.6f\n", SmartIntegral(N, Y, X, x0, 0.5, &f, 300, 1e-8) ); 

  return 0;
}*/

double x0 = 0.9;
double f(double x)
{
  return DOS(DOStypes::SquareLattice, 0.5, x, 0.0)*1.0/(exp(x-x0)-1.0);;
}

int main()
{
  double shift = 0.0;


  GRID grid(300, 5500, 300, 0.9, 0.01, 1e-13, 1.0);
  int N = grid.get_N();

  double* X = new double[N];
  grid.assign_omega(X);
  
  double* Y = new double[N];
  for(int i=0; i<N; i++) Y[i] = f(X[i]);
  PrintFunc("raw",N,Y,X);
  printf("TrapezIntegral: %.6f\n", TrapezIntegral(N, Y, X) ); 

  printf("SmartIntegral: %.6f\n", SmartIntegral(N, Y, X, x0, 0.01, &f, 300, 1e-8) ); 

  return 0;
}

