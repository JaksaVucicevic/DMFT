#include <cstdio>
#include "../source/Input.h"

int main()
{
  Input input("params");
  
  double A = 0;
  double B = 0;
  char C = 'n';
  int D = 0;
  bool E = false;
  int F = 0;
  int N = 4;
  double* a;
  int* b;

  input.ReadParam(A,"Aviator");
  input.ReadParam(B,"Bold Move");
  input.ReadParam(C,"Censored");
  input.ReadParam(D,"Deviation");
  input.ReadParam(E,"Equipment");
  input.ReadParam(F,"Frequency");

  input.ReadParam(N,"Nelements");
  a = new double [N];
  b = new int [N];

  input.ReadArray(N,a,"double_elems");
  input.ReadArray(N,b,"int_elems");
  printf("A: %le B: %le C: %c D: %d E: %s F: %d\n",A,B,C,D,(E)?"true":"false",F);
  printf("a0: %le a1: %le a2: %le a3: %le\n", a[0],a[1],a[2],a[3]);
  printf("b0: %d b1: %d b2: %d b3: %d\n", b[0],b[1],b[2],b[3]);

  delete [] a;
  delete [] b;
  return 0;
}
