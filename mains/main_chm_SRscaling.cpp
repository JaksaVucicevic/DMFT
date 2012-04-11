#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"


void PrintReport(const char* ReportFN, double U, double T, const char* message)
{
  FILE* ReportFile = fopen(ReportFN,"a");
  fprintf(ReportFile,"T=%.3f dU=%.3f %s\n", T, U, message);
  fclose(ReportFile);  
}

int main()
{
  int NdUs = 13;
  double dUmax = 0.3;
  double dUstep = dUmax*2.0/(NdUs-1);

  GRID grid("params");
  
  Result result(&grid);

  grid.assign_omega(result.omega);
  
  double t = 0.5;  
  result.n = 0.5; 
  InitDelta( DOStypes::SemiCircle,	
             grid.get_N(), 		
             0.5, 0.0, 0.01, t, 	
             result.omega, result.Delta); 
    
  InitDOS( DOStypes::SemiCircle, 		
           t, 					
           grid.get_N(), result.omega, result.NIDOS);


  CHM chm("params");

  int N;
  int M;
  double** UT;
  ReadFunc("InstabilityLine", N, M, UT);
  
  printf("N=%d\n",N);

  FILE* SRFile;  
  
  for (int i=0; i<N; i++)
  { 
    double T = UT[1][i];
    double U = UT[0][i];
    
    SRFile = fopen("ScatteringRates","a");  
    fprintf(SRFile,"%.15le",T);
 
    result.PrintResult("initial");

    for (int j  = 0; j<NdUs; j++)
    //for (double dU = -dUmax; dU<dUmax+0.01; dU+=dUstep)
    { 
      double dU = -dUmax + j*dUstep;
      printf("T=%.3f U=%.3f dU=%.3f U+dU=%.3f\n",T, U ,dU, dU+U);
      chm.SetParams(U+dU,T,t);    

      char FN[50];
      sprintf( FN, "CHM.T%.3f.dU%s%.3f", T, (dU>=0) ? "+" : "-", abs(dU) );

      PrintReport("report",  dU,  T, "Working..."); 
      chm.Run(&result);
      PrintReport("report",  dU,  T, "Done!"); 
           
      result.PrintResult(FN);

      fprintf(SRFile," %.15le",imag(result.Sigma[grid.get_N()/2])); 
    }  
    fprintf(SRFile,"\n");
    fclose(SRFile);
  }
  return 0;
}
