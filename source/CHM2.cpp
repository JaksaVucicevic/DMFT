#include <cstdio>
#include "CHM2.h"
#include "GRID.h"
#include "Result.h"
#include "Input.h"
#include "SIAM2.h"
#include "routines.h"
#include <omp.h>

class SIAM2;

void CHM2::Defaults() 
{
    U = 2.0;
    T = 0.05;
    UseBethe = false;
    SIAMeta = 1e-5;
    t = 0.5;
    SiamNt = 8;
    siam = new SIAM2();

    SIAMUseLatticeSpecificG = false;
    LatticeType = DOStypes::SemiCircle;
}

CHM2::CHM2() : Loop()
{
  Defaults();
}

CHM2::CHM2(const char* ParamsFN) : Loop(ParamsFN)
{
  Defaults();
  this->ParamsFN.assign(ParamsFN);
  cout << "-- INFO -- CHM: Params File Name set to:" << this->ParamsFN << endl;

  Input input(ParamsFN);
  
  input.ReadParam(U,"CHM2::U");
  input.ReadParam(T,"CHM2::T");
  input.ReadParam(UseBethe,"CHM2::UseBethe");
  input.ReadParam(SIAMeta,"CHM2::SIAMeta");
  input.ReadParam(t,"CHM2::t");

  input.ReadParam(SIAMUseLatticeSpecificG,"CHM2::SIAMUseLatticeSpecificG");
  input.ReadParam(LatticeType,"CHM2::LatticeType");

  input.ReadParam(SiamNt,"CHM2::SiamNt"); 

  siam = new SIAM2(ParamsFN);
}

void CHM2::ReleaseMemory()
{
  printf("CHM2 release\n");
  siam->~SIAM2();
}

CHM2::~CHM2()
{
  ReleaseMemory();
}

void CHM2::SetParams(double U, double T, double t)
{
  this->U = U;
  this->T = T;
  this->t = t;
}


void CHM2::SetSIAMeta(double SIAMeta)
{
  this->SIAMeta = SIAMeta;
}

void CHM2::SetSIAMUseLatticeSpecificG(bool SIAMUseLatticeSpecificG)
{
  this->SIAMUseLatticeSpecificG = SIAMUseLatticeSpecificG;
}
  
bool CHM2::SolveSIAM()
{
  //SIAM2 siam;
  siam->SetUTepsilon(U,T,0.0);
  if (UseBethe) siam->FORMULA_FOR_G = FormulasForG::NoLattice;
  siam->SetBroadening(SIAMeta);
  //siam->SetUseLatticeSpecificG(SIAMUseLatticeSpecificG, t, LatticeType); 
#ifdef _OMP
  omp_set_num_threads(SiamNt);
#endif
  return siam->Run(r);
}

void CHM2::CalcDelta()
{
  r->PrintResult("CHMDeltaIN");
  
  if (UseBethe) printf("\n\n-- INFO -- CHM2 Self Consistency - Delta = t^G, t=%f\n\n",t);
  if (UseBethe)
    #pragma omp parallel for
    for (int i=0; i<N; i++)        
      r->Delta[i] = sqr(t) * r->G[i];        
  else    
    #pragma omp parallel for
    for (int i=0; i<N; i++) 
      r->Delta[i] = r->omega[i] + r->mu - r->Sigma[i] - 1.0/r->G[i];

  r->PrintResult("CHMDeltaOUT");
  
}
