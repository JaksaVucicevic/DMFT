#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"

double FindLambda(int N, double* lambdas)
{
  int n=5;
  double sigma_min=1e+300;
  double best_lambda;
  for(int i=0; i<N-n; i++)
  { double sum=0;
    for(int j=i; j<i+n; j++)
      sum+=lambdas[j];
    double avg_lambda = sum/n;
    sum=0.0;
    for(int j=i; j<i+n; j++)
      sum+=sqr(lambdas[j]-avg_lambda);
    double sigma=sqrt(sum)/n;
    if (sigma<sigma_min)
    {
      sigma_min=sigma;
      best_lambda = avg_lambda; 
    }
  }
  return best_lambda;
}

double CalcLambda(int N, complex<double>* X0, complex<double>* X1, complex<double>* X2)
{
  double num = 0.0;
  double denom = 0.0;
  for (int i=0; i < N; i++)
  {
    num +=    abs(  ( X0[i] - X1[i] ) 
                   *( X1[i] - X2[i] ) );
    denom += sqr(abs( (X0[i] - X1[i]) ));
  }
  return 1.0 - num/denom;
}

int main(int argc, char* argv[])
{
  /*if (argc!=2) exit(0);
  double T = atof(argv[1]);
  printf(" T=%.3f\n",T);*/
  system("rm lambda.T*");
  
  GRID grid("params");
  int N= grid.get_N();
  Result result1(&grid);
  Result result2(&grid);
  Result result3(&grid);

 
  for (double T=0.30; T>0.05; T-=0.01)
  { for (double U=1.5; U<4.0; U+=0.1)
    { char FN[300];
      sprintf( FN, "CHM.U%.3f.T%.3f", U, T);
      char FN2[300];
      sprintf( FN2, "%s.FAILED", FN);
      if ((not FileExists(FN))and(not FileExists(FN2))) continue;

      char ldFN[300];
      sprintf(ldFN,"lambdas_and_diffs.U%.3f.T%.3f",U,T);
      FILE* ldFile = fopen(ldFN,"w");

      double lambdas[100]; 
      int counter=0;  
      for(int it=1; it<100; it++) 
      { char FN[300];
        sprintf( FN, "CHM.U%.3f.T%.3f.it%d", U, T, it);
        if(not result1.ReadFromFile(FN)) break;
        sprintf( FN, "CHM.U%.3f.T%.3f.it%d", U, T, it+1);
        if(not result2.ReadFromFile(FN)) break;
        sprintf( FN, "CHM.U%.3f.T%.3f.it%d", U, T, it+2);
        if(not result3.ReadFromFile(FN)) break;
      
        double sum = 0;
        for(int i=0; i<N; i++)
          sum +=   sqr( real(result2.G[i]-result1.G[i])) 
                 + sqr( imag(result2.G[i]-result1.G[i]));
        double diff = sqrt(sum)/(2.0*N);
        double simple_diff = abs( imag(result2.G[N/2]-result1.G[N/2]) );
        double lambda = CalcLambda(N, result1.G, result2.G, result3.G);
        fprintf(ldFile,"%d %.15le %.15le %.15le\n",it,diff,lambda, simple_diff);
        lambdas[it-1]=lambda;
        counter++;

        complex<double>* dG = new complex<double>[N];
        for(int i=0; i<N; i++)
          dG[i] = 2.0*((double) N)*(result2.G[i]-result1.G[i])/diff;
        sprintf( FN, "dG.U%.3f.T%.3f.it%d", U, T, it+1);
        PrintFunc(FN,N,dG,result1.omega);
      }      
      fclose(ldFile);
      
      char lFN[300];
      sprintf(lFN,"lambda.T%.3f",T);
      FILE* lFile = fopen(lFN,"a");
      fprintf(lFile,"%.15le %.15le\n",U,FindLambda(counter,lambdas) );
      fclose(lFile);
    }
  } 
  return 0;
}
