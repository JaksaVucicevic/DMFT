#include "../source/Result.h"
#include "../source/GRID.h"
#include "../source/routines.h"
#include "../source/SIAM.h"

int main()
{
  //samo da vidim dal radi git
  //proba za push
  //Milos
  GRID grid("params");

  Result result(&grid);

  InitDelta(DOStypes::SemiCircle, grid.get_N(), 0.1, 0.0, 0.01, 0.5, result.omega, result.Delta);
  InitDOS(DOStypes::SemiCircle, 0.5, grid.get_N(), result.omega, result.NIDOS);
  result.PrintResult("input.res");
  result.mu = 0;
  
  SIAM siam("params");
  
  siam.Run(&result);

  result.PrintResult("output.res");
  
  return 0;
}
