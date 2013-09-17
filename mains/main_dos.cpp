#include "../source/Result.h"
#include "../source/GRID.h"
#include "../source/routines.h"

int main()
{
  GRID grid(0,300,3.5,0.0,0.0);

  Result result(&grid);

  grid.assign_omega(result.omega);

  InitDOS(DOStypes::CubicLattice, 0.5, grid.get_N(), result.omega, result.NIDOS);
  
  PrintFunc("CubicDos",grid.get_N(),result.NIDOS,result.omega);
  
  return 0;
}
