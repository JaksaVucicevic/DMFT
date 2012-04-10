#include "CHM.h"
#include "routines.h"
#include "Broyden.h"


CHM::CHM()
{ 
  SetMaxItsAccrAndCoefs(300, 1e-7, 2, (const int []) {1,1});
  SetPrintIntermediate(false);
  SetHaltOnIterations(false);
  SetBroyden(false, false, 0);
  StartFromInsulator = false;
  StartFromPrevious = false;
  SetUseBethe(false);
  SIAMeta = 5e-5;
  Initialized = false;
}

CHM::~CHM()
{
 if (Initialized) 
 { delete [] LastDelta;
   delete [] omega;
 }
}

void CHM::SetGrid(GRID * grid)
{
  this->grid = grid;
  grid->GetGrid(N,omega);
  if (Initialized) delete [] LastDelta;
  LastDelta = new complex<double>[N];
  Initialized = true;
}

void CHM::SetMaxItsAccrAndCoefs(int MAX_ITS, double Accr, int NtoMix, const int * Coefs)
{
  this->MAX_ITS = MAX_ITS;
  this->Accr = Accr;
  this->Coefs = Coefs;
  this->NtoMix = NtoMix;
}

void CHM::SetPrintIntermediate(bool PrintIntermediate)
{
  this->PrintIntermediate = PrintIntermediate;
}

void CHM::SetHaltOnIterations(bool HaltOnIterations)
{
  this->HaltOnIterations = HaltOnIterations;
}

void CHM::SetBroyden(bool UseBroyden, bool ForceBroyden, double BroydenStartDiff)
{
  this->UseBroyden = UseBroyden;
  this->ForceBroyden = ForceBroyden;
  this->BroydenStartDiff = BroydenStartDiff;
}

void CHM::SetStartFromInsulator(bool StartFromInsulator)
{
 this->StartFromInsulator = StartFromInsulator;
}

void CHM::SetStartFromPrevious(bool StartFromPrevious)
{
 this->StartFromPrevious = StartFromPrevious;
}

void CHM::SetUseBethe(bool UseBethe)
{
 this->UseBethe = UseBethe;
}

void CHM::SetSIAMeta(double eta)
{
 this->SIAMeta = eta;
}
 


void CHM::Run(double n, double U, double T, double t, double &mu, int DOStype, const char* FNDos)
{
  if (!Initialized) { printf("CHM NOT initialized!"); exit(1); }
   
  //relevant quantities
  complex<double>* Delta = new complex<double>[N];	//hybridization function
  complex<double>* G = new complex<double>[N];		//greens function
  complex<double>* Sigma = new complex<double>[N];	//self-energy
  //complex<double>* dos = new complex<double>[N];  


  //Initialize mixer
  Mixer< complex<double> > mixer(N, NtoMix, Coefs, (UseBroyden) ? BroydenStartDiff : Accr);


  //----- Delta initialization ------//
  for(int i=0; i<N; i++) 
               Delta[i] = 2.0 * sqr(t) /
                         (    omega[i] 
                            + ((omega[i] < -2.0 * t) ? -1.0 : +1.0) * sqrt( complex<double>( sqr(omega[i]) - 4.0 * sqr(t), 0.0) )
                         ) ;
  PrintFunc("DeltaBre", N, Delta, omega);
/*  if (StartFromPrevious)
    for(int i=0; i<N; i++) Delta[i] = LastDelta[i];
  else 
  {
    if (StartFromInsulator)
      //InitDelta(DOStypes::Insulator, N,  1.0, 0.0, 0.02, t, omega, Delta);
      InitInsulatorDelta(U,N,t,omega,Delta);
    else
    {  //InitDelta(DOStype, N,  1.0, 0.0, 0.02, t, omega, Delta, FNDos);
       if (UseBethe)
          for(int i=0; i<N; i++) 2.0 * sqr(t) /
               (    omega[i] 
                  - sqrt( complex<double>( sqr(omega[i]) + 4.0 * sqr(t), 0.0) )
               ) ;
           
    }
  }
*/

  mixer.Mix(Delta);
  //---------------------------------//

  //---------init SIAM---------------//
  SIAM siam;
  siam.Initialize(grid);  
  siam.InitImpurity(U, T, 0);
  siam.SetDOStype_CHM(DOStype, t, FNDos);
  siam.SetBroadening(SIAMeta);
  if (DOStype==DOStypes::SemiCircle and UseBethe) 
  {  siam.SetIsBethe(true);
     printf("NOTE: Using Bethe simplified Delta = t^2 G\n");
  }
  //---------------------------------//
 
  //initialize broyden
  Broyden B;
  B.SetParameters(N, MAX_ITS, 1.0, 0.01, Accr);

  // Broyden status: 0 - Waiting for mixer to reach BroydenStartDiff
  //                 1 - Running
  //                 2 - Suspended
  int BroydenStatus = 0;

  //Halt on first iteration if HaltOnIterations
  int Halt = (HaltOnIterations) ? 1 : 0; 

  printf("*******************************************************************\n");
  printf(" DMFT loop: Clean Hubbard Model, n=%.3f U = %.3f, T = %.3f\n", n, U, T);
  printf("*******************************************************************\n");

  //------------ DMFT loop-------------//
  for (int it = 1; it<=MAX_ITS; it++)
  {  printf("--- DMFT Iteration %d ---\n", it);
   
    //set accr for siam broyden
    /* siam.SetBroydenParameters(100, (BroydenStatus == 1) ? max(B.CurrentDiff, 1e-6) 
                                                         : min(mixer.CurrentDiff,BroydenStartDiff) );
    */
    
    printf("Integral Delta = %.6f\n",imag(TrapezIntegral(N,Delta,omega)));
    
     //----- solve SIAM ------//
     if ( siam.Run_CHM(n, Delta,  G, Sigma, mu)
          and (BroydenStatus == 1) and (!ForceBroyden) ) //if clipping, turn off broyden
     {   BroydenStatus++; 
         mixer.Initialize(N, NtoMix, Coefs, Accr);
         mixer.Mix(Delta);
     } 
     //-----------------------//

     //halt
     if (it==Halt)
     {
       siam.PrintResults("intermediate");
       printf("Next stop: ");
       cin >> Halt; 
     }

     //print out intermediate results
     if (PrintIntermediate)
     {  char FN[50];
        sprintf(FN,"CHM.n%.3f.U%.3f.T%.3f.it%d", n, U, T, it);
        siam.PrintResults(FN);
     }

     //--- self-consistency ---// 
     if (DOStype==DOStypes::SemiCircle and UseBethe)
     {  for (int i=0; i<N; i++)        
          Delta[i] = sqr(t)*G[i];        
        printf("Delta = t^2 G\n");
     }
     else    
       for (int i=0; i<N; i++) 
       {  Delta[i] = omega[i] + mu - Sigma[i] - complex<double>(1.0)/G[i];
          complex<double> d = omega[i] + mu - Sigma[i] - conj(G[i])/(sqr(real(G[i]))+sqr(imag(G[i])));
/*          if ((i>N/2-50)and(i<N/2+50)) printf("omega: %.5f mu: %.5f Sigma: (%.5f, %.5f) G: (%.5f, %.5f) 1/G:( %.5f, %.5f) Delta: (%.5f, %.5f) (%.5f, %.5f)\n",
                                          omega[i], mu, real(Sigma[i]), imag(Sigma[i]), real(G[i]), imag(G[i]), real(complex<double>(1.0)/G[i]),
                                                            imag(complex<double>(1.0)/G[i]), real(Delta[i]), imag(Delta[i]), real(d), imag(d) ); */
       }
         
    //------------------------//

     // now mix and check if converged
     bool conv = false;
     if (BroydenStatus == 1) 
       conv = B.CalculateNew(Delta,it);
     else
     { if (mixer.Mix(Delta))
         if ((UseBroyden)and(BroydenStatus == 0)) 
         { B.TurnOn(it); //switch to broyden if mixer converged
           BroydenStatus++;
           B.CurrentDiff = mixer.CurrentDiff;
         } 
         else conv = true;
     }
     for (int i=0; i<N; i++) LastDelta[i] = Delta[i];
     if (conv) { IterationsMade = it; break; }
  }
  //-----------------------------------//
  
  if (BroydenStatus == 1) B.TurnOff();

  //print out final results
  char FN[50]; 
  sprintf(FN,"CHM.N%d.eta%le.n%.3f.U%.3f.T%.6f", N, SIAMeta, n, U, T);
  siam.PrintResults(FN);
  //siam.PrintModel();

  //release memory
  delete [] Delta;
  delete [] G;
  delete [] Sigma;
}

void CHM::InitDeltaFromFile(const char* FNDelta, int Nlog, int Nlin, double omega_lin_max, double omega_max, double omega_min)
{ //----------------- OBRNUTI N i M !!!!!! --------------------//
  GRID grid;
  grid.InitGrid(Nlog,Nlin,omega_lin_max,omega_max,omega_min);
  double* omega;
  int N, M;
  grid.GetGrid(N,omega);

  double** d;
  ReadFunc(FNDelta, N, M, d);
  complex<double>* Delta = new complex<double>[N];
  for (int i=0; i<N; i++) Delta[i] = complex<double>(d[i][2], d[i][3]);
  
  for (int i=0; i<this->N; i++)
    LastDelta[i] = grid.interpl(Delta,this->omega[i]);
  
  for  (int i=0; i<N; i++)
    delete [] d[i];
  delete [] d;
  delete [] Delta;
  
  StartFromPrevious = true;
}
                    
