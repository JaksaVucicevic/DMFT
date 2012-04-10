#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"

void PrintReport(const char* ReportFN, double U, double T, const char* message)
{
  FILE* ReportFile = fopen(ReportFN,"a");
  fprintf(ReportFile,"T=%.3f U=%.3f %s\n", T, U, message);
  fclose(ReportFile);  
}

int main()
{
  double Ustart = 2.5;
  double Uend = 4.0;
  double Ustep_max = 0.5;
  double Ustep_min = 0.0125;

  double Tstart;
  double Tend;
  double Tstep;

  Input input("params");
  input.ReadParam(Tstart, "main::Tstart");
  input.ReadParam(Tend, "main::Tend");
  input.ReadParam(Tstep, "main::Tstep");
  GRID grid("params");
  
  Result result(&grid);

  grid.assign_omega(result.omega);
  
  double t = 0.5;
  input.ReadParam(t,"CHM::t"); 
  
  InitDOS( DOStypes::SemiCircle, 		
           t, 					
           grid.get_N(), result.omega, result.NIDOS);

  double n = 0.5;
  input.ReadParam(n,"main::n");
  result.n = n; 

  CHM chm("params");
  double Ustep;

  for (double T=Tstart; T<Tend; T+=Tstep)
  { 
    InitDelta( DOStypes::SemiCircle,	
               grid.get_N(), 		
               0.5, 0.0, 0.01, t, 	
               result.omega, result.Delta); 
    
    InitDOS( DOStypes::SemiCircle, 		
             t, 					
             grid.get_N(), result.omega, result.DOS);

    double LastDOS0 = result.DOS[grid.get_N()/2];      

    double Ustep = Ustep_max; 
    double U = Ustart;
    
    bool ReachedEnd = false;
    do 
    { 
      chm.SetParams(U,T,t);
      
      Result resCopy(result);

      char FN[50];
      sprintf( FN, "CHM.%s.T%.3f.U%.3f", (Ustep>0.0) ? "FromMet" : "FromIns", T, U );

      FILE* f = fopen(FN,"r");
      if (f!=NULL) 
        result.ReadFromFile(FN);
      else
      {  
         PrintReport("report",  U,  T, "Working..."); 
         chm.Run(&result);
         PrintReport("report",  U,  T, "Done!"); 
      }
      
      result.PrintResult(FN);
  
      if ( ( (Ustep > 0.0) and (LastDOS0 / result.DOS[grid.get_N()/2] > 3.0) )
         )
      { 
           PrintReport("report",  U,  T, "Insultor Found");            

           result.CopyFrom(resCopy);
           U -= Ustep;  
           Ustep *= 0.5;

           if (abs(Ustep)<Ustep_min) 
           {  
              Ustep = - Ustep_max; 
              U = Uend + Ustep_max;
              
              InitDelta( DOStypes::SemiCircle,	
                         grid.get_N(), 		
                         0.5, 0.0, 0.01, t, 	
                         result.omega, result.Delta);
              InitDOS( DOStypes::SemiCircle, 		
                       t, 					
                       grid.get_N(), result.omega, result.DOS);
           }
           
      }
      else
      if ( ( (Ustep < 0.0) and (result.DOS[grid.get_N()/2] / LastDOS0 > 3.0) ) 
         )
      { 
           PrintReport("report",  U,  T, "Metal Found");            

           result.CopyFrom(resCopy);
           U -= Ustep;  
           Ustep *= 0.5;
                   

           if (abs(Ustep)<Ustep_min) ReachedEnd = true;
      }

      LastDOS0 = result.DOS[grid.get_N()/2];      

      U+=Ustep;

    } while (!ReachedEnd);  
  }
  return 0;
}
