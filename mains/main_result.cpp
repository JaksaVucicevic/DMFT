#include "../source/Result.h"
#include "../source/GRID.h"
#include "../source/routines.h"

int main()
{
  // initialize grid. params are read from input file
  GRID grid(0,500,2.0,0.0,0.0);

  // initialize a result object on grid
  Result result(&grid);

  // fill in Delta by initializing it to semi-circle
  InitDOS(DOStypes::Cubic, 0.5, grid.get_N(), result.omega, result.NIDOS)

  PrintFunc("CubicDos",result.NI
  result.mu = 0.1;
  result.n = 0.5;
  result.mu0 = 0.3; 
 
  //make a copy
  Result r2(result);
  //r2.CopyFrom(&result);
  r2.mu = 0.2;  
  r2.Delta[1000] = 7;

  result.PrintResult("Result.dat");
  r2.PrintResult("Result.dat.c");  

  Result r3(&grid);
  // this way mu, n and mu0 are not read
  r3.ReadFromFile("Result.dat.c");
  r3.n = 0.3;  
  r3.PrintResult("Result.dat.cc");  
  
  return 0;
}
