#include "../source/Result.h"
#include "../source/GRID.h"
#include "../source/routines.h"
#include "../source/SIAM.h"

int main()
{
  //samo da vidim dal radi git
  //proba za push
  
  GRID grid("params");

  Result result(&grid);

  InitDelta(DOStypes::SemiCircle, grid.get_N(), 0.1, 0.0, 0.01, 0.5, result.omega, result.Delta);
  result.mu = 0;
  
  SIAM siam("params");
  
  siam.Run(&result);

  result.PrintResult("Result.dat");
  
  return 0;
}
