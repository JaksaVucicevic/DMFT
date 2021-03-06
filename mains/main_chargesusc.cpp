#include "../source/CHM.h"
#include "../source/GRID.h"
#include "../source/Result.h"
#include "../source/routines.h"
#include "../source/Input.h"

/*
int main()
{

  GRID grid("params");

  Result result(&grid);

  grid.assign_omega(result.omega);

 
  for (double U=0.05; U<4.45; U+=0.1)
  {

    for (double T=0.002; T<0.12; T+=0.01)
    {
      printf("U,T: %.3f, %.3f\n",U,T);
      char ResFN[100]; 
      sprintf( ResFN, "/nfs/jaksa/WL_IPT/run/CHM.U%.3f.T%.3f%s", U, T, ((U>=2.498)and(U<=3.5)) ? ".FromSC" : "");
    
      bool b = result.ReadFromFile(ResFN);
      if (!b) 
      { 
         printf("FILE NOT FOUND!\n"); 
         continue;
      }

      double chi1=0, chi3=0;
      result.ChargeSusceptibility(T,chi1,chi3);

      FILE* chargeSuscFile = fopen("/home/jaksa/WL_ChargeSusc/NEW3/chargeSuscUT","a");
      fprintf(chargeSuscFile,"%.15le %.15le %.15le %.15le\n", U, T, chi1, chi3);
      fclose(chargeSuscFile);

      char csTFN[100];
      sprintf( csTFN, "/home/jaksa/WL_ChargeSusc/NEW3/chargeSusc.T%.3f", T);
      FILE* csTFile = fopen(csTFN,"a");
      fprintf(csTFile,"%.15le %.15le %.15le\n", U, chi1, chi3);
      fclose(csTFile);

      char csUFN[100];
      sprintf( csUFN, "/home/jaksa/WL_ChargeSusc/NEW3/chargeSusc.U%.3f", U);
      FILE* csUFile = fopen(csUFN,"a");
      fprintf(csUFile,"%.15le %.15le %.15le\n", T, chi1, chi3);
      fclose(csUFile);
    }
    FILE* chargeSuscFile = fopen("/home/jaksa/WL_ChargeSusc/NEW3/chargeSuscUT","a");
    fprintf(chargeSuscFile,"\n");
    fclose(chargeSuscFile);
  }
  
  return 0;
}

*/




/*
int main()
{

  GRID grid("params");

  Result result(&grid);

  grid.assign_omega(result.omega);

  Input input("params");
  double t = 0.5; // always set a default value in case parameter is not found in the input file!!!
  input.ReadParam(t,"CHM::t");
   
  InitDOS( DOStypes::SemiCircle, 			// type of DOS 
           t, 					// hopping amplitude
           grid.get_N(), result.omega, result.NIDOS);	// total number of omega grid points, omega-grid array, the array in which NIDOS will be stored


  CHM chm("params");

  for (double U=0.1; U<3.5; U+=0.3)
  {
    InitDelta( DOStypes::SemiCircle,	// type of DOS 
             grid.get_N(), 		// total number of omega grid points 
             0.5, 0.0, 0.01, t, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping amplitude
             result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored
     
    for (double T=0.2; T>0.01; T-=0.03)
    {
      result.n = 0.50;
      result.mu0 = 0.0;
  
      chm.SetParams(U,T,t);
      chm.Run(&result);

      double chi1=0, chi3=0;
      result.ChargeSusceptibility(T,chi1,chi3);

      FILE* chargeSuscFile = fopen("chargeSuscUT","a");
      fprintf(chargeSuscFile,"%.15le %.15le %.15le %.15le\n", U, T, chi1, chi3);
      fclose(chargeSuscFile);

      char csTFN[50];
      sprintf( csTFN, "chargeSusc.T%.3f", chm.get_T());
      FILE* csTFile = fopen(csTFN,"a");
      fprintf(csTFile,"%.15le %.15le %.15le\n", U, chi1, chi3);
      fclose(csTFile);

      char csUFN[50];
      sprintf( csUFN, "chargeSusc.U%.3f", chm.get_U());
      FILE* csUFile = fopen(csUFN,"a");
      fprintf(csUFile,"%.15le %.15le %.15le\n", T, chi1, chi3);
      fclose(csUFile);
    }
    FILE* chargeSuscFile = fopen("chargeSuscUT","a");
    fprintf(chargeSuscFile,"\n");
    fclose(chargeSuscFile);
  }
  
  return 0;
}

*/















int main()
{

  GRID grid("params");

  Result result(&grid);

  grid.assign_omega(result.omega);

  Input input("params");
  double t = 0.5; // always set a default value in case parameter is not found in the input file!!!
  input.ReadParam(t,"CHM::t");
   
  InitDOS( DOStypes::SemiCircle, 			// type of DOS 
           t, 					// hopping amplitude
           grid.get_N(), result.omega, result.NIDOS);	// total number of omega grid points, omega-grid array, the array in which NIDOS will be stored


  CHM chm("params");

  for (double U=0.1; U<4.0; U+=0.05)
  {
    InitDelta( DOStypes::SemiCircle,	// type of DOS 
             grid.get_N(), 		// total number of omega grid points 
             0.5, 0.0, 0.01, t, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping amplitude
             result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored
     
    for (double T=0.15; T>=0.009; T-=0.01)
    {
      result.n = 0.52;
      result.mu0 = 0.0;
  
      chm.SetParams(U,T,t);
      chm.SetParams(U,T,t);
      bool notconverged = chm.Run(&result);
  
      char FN[50];
      sprintf( FN, "CHM.U%.3f.T%.3f.n%.3f%s", chm.get_U(), chm.get_T(), result.n, (notconverged) ? ".FAILED" : "");
      result.PrintResult(FN);
  
      FILE* chargeSuscFile = fopen("chargeSuscUT","a");
      fprintf(chargeSuscFile,"%.15le %.15le %.15le %.15le\n", U, T, 0.02/(result.mu - U/2.0),result.mu);
      fclose(chargeSuscFile);

      char csTFN[50];
      sprintf( csTFN, "chargeSusc.T%.3f", chm.get_T());
      FILE* csTFile = fopen(csTFN,"a");
      fprintf(csTFile,"%.15le %.15le %.15le\n", U, 0.02/(result.mu - U/2.0), result.mu);
      fclose(csTFile);

      char csUFN[50];
      sprintf( csUFN, "chargeSusc.U%.3f", chm.get_U());
      FILE* csUFile = fopen(csUFN,"a");
      fprintf(csUFile,"%.15le %.15le %.15le\n", T, 0.02/(result.mu - U/2.0), result.mu);
      fclose(csUFile);
    }
    FILE* chargeSuscFile = fopen("chargeSuscUT","a");
    fprintf(chargeSuscFile,"\n");
    fclose(chargeSuscFile);
  }
  
  return 0;
}


/*
int main()
{

  GRID grid("params");

  Result result(&grid);

  grid.assign_omega(result.omega);

  Input input("params");
  double t = 0.5; // always set a default value in case parameter is not found in the input file!!!
  input.ReadParam(t,"CHM::t");
  

 

  InitDOS( DOStypes::SemiCircle, 			// type of DOS 
           t, 					// hopping amplitude
           grid.get_N(), result.omega, result.NIDOS);	// total number of omega grid points, omega-grid array, the array in which NIDOS will be stored


  CHM chm("params");


  InitDelta( DOStypes::SemiCircle,	// type of DOS 
               grid.get_N(), 		// total number of omega grid points 
               0.5, 0.0, 0.01, t, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping amplitude
               result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored

  for (double U=4.5; U<4.51; U+=0.05)
  for (double T=0.06; T<0.08; T+=0.01)
  for (double n=0.5; n<0.8; n+=0.025)
  {
   
    result.n = n;
    result.mu0 = 0;
    chm.SetParams(U,T,t);
    bool notconverged = chm.Run(&result);
  
    char FN[50];
    sprintf( FN, "CHM.U%.3f.T%.3f.n%.3f%s", chm.get_U(), chm.get_T(), n, (notconverged) ? ".FAILED" : "");
    result.PrintResult(FN);

    char mu_vs_nFN[50];
    sprintf( mu_vs_nFN, "mu_vs_n.U%.3f.T%.3f", chm.get_U(), chm.get_T());
    FILE* mu_vs_nFile = fopen(mu_vs_nFN,"a");
    fprintf(mu_vs_nFile,"%.15le %.15le\n", result.n, result.mu);
    fclose(mu_vs_nFile);

    if (notconverged)
        InitDelta( DOStypes::SemiCircle,	// type of DOS 
               grid.get_N(), 		// total number of omega grid points 
               0.5, 0.0, 0.01, t, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping amplitude
               result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored

  }

  
  return 0;
}
*/

/*
int main()
{

  GRID grid("params");

  Result result(&grid);

  grid.assign_omega(result.omega);

  Input input("params");
  double t = 0.5; // always set a default value in case parameter is not found in the input file!!!
  input.ReadParam(t,"CHM::t");
  

 

  InitDOS( DOStypes::SemiCircle, 			// type of DOS 
           t, 					// hopping amplitude
           grid.get_N(), result.omega, result.NIDOS);	// total number of omega grid points, omega-grid array, the array in which NIDOS will be stored


  CHM chm("params");

  double T=0.01;
  double n = 0.52;
  for (double U=0.1; U<3.5; U+=0.3)
  {  
    InitDelta( DOStypes::SemiCircle,	// type of DOS 
               grid.get_N(), 		// total number of omega grid points 
               0.5, 0.0, 0.01, t, 	// hybridization-V, chemical potential-mu (defines the omega shift), broadening eta, and hopping amplitude
               result.omega, result.Delta);  // omega-grid array and the array in which Delta will be stored

    result.n = n;
    result.mu0 = 0;
    chm.SetParams(U,T,t);
    bool notconverged = chm.Run(&result);
  
    char FN[50];
    sprintf( FN, "CHM.U%.3f.T%.3f.n%.3f%s", chm.get_U(), chm.get_T(), n, (notconverged) ? ".FAILED" : "");
    result.PrintResult(FN);

    char csUFN[50];
    sprintf( csUFN, "chargeSusc.T%.3f", chm.get_T());
    FILE* csUFile = fopen(csUFN,"a");
    fprintf(csUFile,"%.15le %.15le %.15le\n", U, 0.02/(result.mu - U/2.0), result.mu);
    fclose(csUFile);

  }

  
  return 0;
}
*/
