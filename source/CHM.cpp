#include <cstdio>
#include "CHM.h"
#include "GRID.h"
#include "Result.h"
#include "Input.h"
#include "SIAM.h"
#include "routines.h"
#include <omp.h>

class SIAM;

void CHM::Defaults() 
{
    U = 2.0;
    T = 0.05;
    UseBethe = false;
    SIAMeta = 1e-5;
    UseSmartSIAMeta = false;
    UseFixedMuSIAMRun = false;

    t = 0.5;
    SiamNt = 8;
    siam = new SIAM();

    SIAMUseLatticeSpecificG = false;
    LatticeType = DOStypes::SemiCircle;
}

CHM::CHM() : Loop()
{
  Defaults();
}

CHM::CHM(const char* ParamsFN) : Loop(ParamsFN)
{
  Defaults();
  this->ParamsFN.assign(ParamsFN);
  cout << "-- INFO -- CHM: Params File Name set to:" << this->ParamsFN << endl;

  Input input(ParamsFN);
  
  input.ReadParam(U,"CHM::U");
  input.ReadParam(T,"CHM::T");
  input.ReadParam(UseBethe,"CHM::UseBethe");
  input.ReadParam(SIAMeta,"CHM::SIAMeta");
  input.ReadParam(UseSmartSIAMeta,"CHM::UseSmartSIAMeta");
  input.ReadParam(UseFixedMuSIAMRun,"CHM::UseFixedMuSIAMRun");

  input.ReadParam(t,"CHM::t");

  input.ReadParam(SIAMUseLatticeSpecificG,"CHM::SIAMUseLatticeSpecificG");
  input.ReadParam(LatticeType,"CHM::LatticeType");

  input.ReadParam(SiamNt,"CHM::SiamNt"); 

  siam = new SIAM(ParamsFN);
}

void CHM::ReleaseMemory()
{
  printf("CHM release\n");
  siam->~SIAM();
}

CHM::~CHM()
{
  ReleaseMemory();
}

void CHM::SetParams(double U, double T, double t)
{
  this->U = U;
  this->T = T;
  this->t = t;
}

void CHM::SetUseBethe(bool UseBethe)
{
  this->UseBethe = UseBethe;
}

void CHM::SetSIAMeta(double SIAMeta, bool UseSmartSIAMeta)
{
  this->SIAMeta = SIAMeta;
  this->UseSmartSIAMeta=UseSmartSIAMeta;
}

void CHM::SetSIAMUseLatticeSpecificG(bool SIAMUseLatticeSpecificG)
{
  this->SIAMUseLatticeSpecificG = SIAMUseLatticeSpecificG;
}
  
bool CHM::SolveSIAM()
{
  //SIAM siam;
  siam->SetUTepsilon(U,T,0.0);
  siam->SetIsBethe(UseBethe);

  double eta = SIAMeta;
  if (UseSmartSIAMeta)
  { if (Iteration>10)
      eta *= 0.3;
    if (Iteration>20)
      eta *= 0.05;
  }
  siam->SetBroadening(eta);

  siam->SetUseLatticeSpecificG(SIAMUseLatticeSpecificG, t, LatticeType); 
#ifdef _OMP
  omp_set_num_threads(SiamNt);
#endif

  if (UseFixedMuSIAMRun)
    return siam->Run(r);
  else
    return siam->Run_CHM(r);
}

void CHM::CalcDelta()
{
  r->PrintResult("CHMDeltaIN");
  
  if (UseBethe) printf("\n\n-- INFO -- CHM Self Consistency - Delta = t^G, t=%f\n\n",t);
  if (UseBethe)
    #pragma omp parallel for
    for (int i=0; i<N; i++)        
      r->Delta[i] = sqr(t) * r->G[i];        
  else    
    #pragma omp parallel for
    for (int i=0; i<N; i++) 
    {  r->Delta[i] = r->omega[i] + r->mu - r->Sigma[i] - 1.0/r->G[i];
       //if (imag(r->Delta[i])>0)  r->Delta[i]=r->omega[i] + r->mu - real(r->Sigma[i]) - real(1.0/r->G[i]);
    }
  r->PrintResult("CHMDeltaOUT");

  if (PrintIntermediate)
  { char FN[300];
    sprintf( FN, "CHM.U%.3f.T%.3f.it%d", U, T, Iteration);
    r->PrintResult(FN);
  }  
  
}
