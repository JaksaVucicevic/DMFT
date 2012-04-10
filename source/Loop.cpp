#include "Broyden.h"
#include "Mixer.h"
#include "Result.h"
#include "Input.h"
#include "GRID.h"
#include "Loop.h"

void Loop::Defaults()
{
    printf("\n\n\n\nLOOP DEFAULTS <<<<<<<<<<<<<<<<<<<<\n\n");
    UseBroyden = true;
    ForceBroyden = false;
    BroydenStartDiff =5e-3;

    //---- Mixer Options ------//
    NtoMix = 2;
    Coefs = new int[NtoMix];
    Coefs[0] = 1;
    Coefs[1] = 0;
    
    //---- Loop Options -------//
    MAX_ITS = 300;
    Accr = 5e-5;

    //---- PrintOut/Debugging optins----//
    PrintIntermediate = false;
    HaltOnIterations = false;
    ForceSymmetry = false;
}

Loop::Loop()
{
  Defaults();
}

Loop::Loop(const char* ParamsFN)
{
  Defaults();
  
  this->ParamsFN.assign(ParamsFN);
  cout << "-- INFO -- Loop: Params File Name set to:" << this->ParamsFN << endl;

  Input input(ParamsFN);
  
  input.ReadParam(UseBroyden,"Loop::UseBroyden");
  input.ReadParam(ForceBroyden,"Loop::ForceBroyden");
  input.ReadParam(BroydenStartDiff,"Loop::BroydenStartDiff");
  input.ReadParam(NtoMix,"Loop::NtoMix");
  delete [] Coefs;
  Coefs = new int[NtoMix];
  Coefs[0] = 1;
  printf("-------- LOOP:: C0 = %d, C1 = %d\n",Coefs[0],Coefs[1]);   
  input.ReadArray(NtoMix, Coefs, "Loop::Coefs");
  input.ReadParam(MAX_ITS,"Loop::MAX_ITS");
  input.ReadParam(Accr,"Loop::Accr");
  input.ReadParam(PrintIntermediate,"Loop::PrintIntermediate");
  input.ReadParam(HaltOnIterations,"Loop::HaltOnIterations");
  input.ReadParam(ForceSymmetry,"Loop::ForceSymmetry");
  
  printf("HOI: %s PI: %s\n", (HaltOnIterations) ? "yes" : "no", (PrintIntermediate) ? "yes" : "no");
}

void Loop::ReleaseMemory()
{
  printf("Loop release\n");
  delete [] Coefs;
}

Loop::~Loop()
{
  ReleaseMemory();
}

//------- Option Stetters --------//

void Loop::SetGrid(GRID* grid)
{
  this->grid = grid;
  N = grid->get_N();
}

void Loop::SetMixerOptions(int NtoMix, const int * Coefs)
{

  delete [] Coefs;
  this->NtoMix = NtoMix;
  this->Coefs = new int[NtoMix];
  for (int i=0; i<NtoMix; i++) this->Coefs[i] = Coefs[i];;
}

void Loop::SetBroydenOptions(bool UseBroyden, bool ForceBroyden, double BroydenStartDiff)
{
  this->UseBroyden = UseBroyden;
  this->ForceBroyden = ForceBroyden;
  this->BroydenStartDiff = BroydenStartDiff;
}

void Loop::SetLoopOptions(int MAX_ITS, double Accr)
{
  this->MAX_ITS = MAX_ITS;
  this->Accr = Accr;
}

void Loop::SetPrintOutOptions(bool PrintIntermediate, bool HaltOnIterations)
{
  this->PrintIntermediate = PrintIntermediate;
  this->HaltOnIterations = HaltOnIterations;
  if (HaltOnIterations) printf("-- INFO -- Loop: Halt on iterations set ON\n");
}

//---------------------------------------------------//

bool Loop::Run(Result* r)
{
  this->r = r;
  this->grid = r->grid;
  N = r->grid->get_N();
  
  //Initialize mixer
  printf("|||||||||||||||||||||||||| LOOP:: C0 = %d, C1 = %d\n",Coefs[0],Coefs[1]);
  Mixer< complex<double> > mixer(N, NtoMix, Coefs, (UseBroyden) ? BroydenStartDiff : Accr);
  mixer.Mix(r->Delta);

  //initialize broyden
  Broyden B;
  B.SetParameters(N, MAX_ITS, 1.0, 0.01, Accr);
  int BroydenStatus = 0;
  // Broyden status: 0 - Waiting for mixer to reach BroydenStartDiff
  //                 1 - Running
  //                 2 - Suspended
  

  //Halt on first iteration if HaltOnIterations
  int Halt = (HaltOnIterations) ? 1 : 0; 

  bool converged = false;
  //------------ DMFT loop-------------//
  for (int it = 1; it<=MAX_ITS; it++)
  {  printf("--- DMFT Loop Iteration %d ---\n", it);
     
    //set accr for siam broyden
    /* siam.SetBroydenParameters(100, (BroydenStatus == 1) ? max(B.CurrentDiff, 1e-6) 
                                                         : min(mixer.CurrentDiff,BroydenStartDiff) );
    */
    
     //----- solve SIAM ------//
     if ( SolveSIAM() ) return true;

// =========     TODO    ========= Handling errors in SIAM
/*     if ( SolveSIAM() and (BroydenStatus == 1) and (!ForceBroyden) ) //if clipping, turn off broyden, unless broyden is forced
     {   BroydenStatus++; 
         mixer.Initialize(N, NtoMix, Coefs, Accr);
         mixer.Mix(r->Delta);
     } */
     //-----------------------//

     //halt
     if (it==Halt)
     {
       r->PrintResult("intermediate");
       printf("Next stop: ");
       cin >> Halt; 
     }

     //print out intermediate results
     if (PrintIntermediate)
     {  char FN[50];
        sprintf(FN,"intermediate.%d", it);
        //sprintf(FN,"intermediate");
        r->PrintResult(FN);
     }

     //--- self-consistency ---// 
     CalcDelta(); 
     //------------------------//

     // check for nans
     //#pragma omp parallel for
     for (int i = 0; i < N; i++) if ( r->Delta[i] != r->Delta[i] ) { printf("nan in Delta!!!!\n"); return true; }

     // clip off
     bool ClippingDelta= false;

     #pragma omp parallel for
     for (int i = 0; i < N; i++) 
       if ( imag(r->Delta[i]) > 0.0 ) 
       {  ClippingDelta = true;
          r->Delta[i] = complex<double>( real(r->Delta[i]), -imag(r->Delta[i]) ); 
       }
     if (ClippingDelta) printf("||||||||||||| Clipping Delta!!!!\n");  
     // force symmetry
     if (ForceSymmetry)
       #pragma omp parallel for 
       for (int i = 0; i < N/2 - 1; i++) r->Delta[i] = -conj(r->Delta[N - 1 - i]);
     


     // now mix and check if converged
     int conv = 0;
     if (BroydenStatus == 1) 
       conv = B.CalculateNew(r->Delta,it);
     else
     { if (mixer.Mix(r->Delta))
         if ((UseBroyden)and(BroydenStatus == 0)) 
         { B.TurnOn(it); //switch to broyden if mixer converged
           BroydenStatus++;
           B.CurrentDiff = mixer.CurrentDiff;
         } 
         else conv = 1;
     }
     if (conv==1) { converged = true; break; }
  }
  //-----------------------------------//
  
  if (BroydenStatus == 1) B.TurnOff();
  return !converged;
}

bool Loop::SolveSIAM()
{
  printf("SS Loop"); 
  return false;
}

void Loop::CalcDelta()
{
  printf("CD Loop");
}
