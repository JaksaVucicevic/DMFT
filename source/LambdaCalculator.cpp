#include "routines.h"
#include "Input.h"
#include "LambdaCalculator.h"
#include <cmath>
#include <cstdio>


LambdaCalculator::LambdaCalculator()
{
  Defaults();

}
LambdaCalculator::LambdaCalculator(int N, int M, double T)
{
  Defaults();
  SetNMT(N,M,T);
}

LambdaCalculator::LambdaCalculator(const char* paramsFN)
{
  Defaults();

  Input input(paramsFN);

  input.ReadParam(N,"LambdaCalculator::N");
  input.ReadParam(M,"LambdaCalculator::M");
  input.ReadParam(T,"LambdaCalculator::T");
  SetNMT(N, M, T);
  offset = N/2;		
  input.ReadParam(offset,           "LambdaCalculator::offset");
  input.ReadParam(Nlambdas,         "LambdaCalculator::Nlambdas");
  input.ReadArray(Nlambdas, Nfreqs, "LambdaCalculator::Nfreqs");
  input.ReadParam(Ndiffs,           "LambdaCalculator::Ndiffs");
  input.ReadArray(Ndiffs, is,       "LambdaCalculator::is");

  input.ReadParam(continued_Nlambdas,                   "LambdaCalculator::continued_Nlambdas");
  input.ReadArray(continued_Nlambdas, continued_Nfreqs, "LambdaCalculator::continued_Nfreqs");
  input.ReadParam(continued_Ndiffs,                     "LambdaCalculator::continued_Ndiffs");
  input.ReadArray(continued_Ndiffs, continued_is,       "LambdaCalculator::continued_is");

  input.ReadParam( returned_Nfreq,         "LambdaCalculator::returned_Nfreq");
  input.ReadParam( returned_FromContinued, "LambdaCalculator::returned_FromContinued");

  input.ReadParam( DoContinued, "LambdaCalculator::DoContinued");
  char* FN = new char[200];
  int err = input.ReadParam( FN, "LambdaCalculator::OutputFN");
  if (err!=-1)  sprintf(OutputFN,"%s",FN);				// output options
  delete [] FN;

  input.ReadParam( DoOutput, "LambdaCalculator::DoOutput");
  input.ReadParam( DoPrintOutContinuedX, "LambdaCalculator::DoPrintOutContinuedX");
}
LambdaCalculator::~LambdaCalculator()
{
  ReleaseMemory();
}

void LambdaCalculator::Defaults()
{ 
  printf("-- info -- LambdaCalculator: Defaults\n");
  Initialized = false;
  counter = 0;


  SetNMT(2000, 1000, 0.05);
  offset = N/2;					//used only for the input X. not used with continued X				
  
  Nlambdas = 2;				//number of different lambdas to be calculated (max 10)
  Nfreqs[0] = 1;				//the number of freqs used for calculation of eachj lambda
  Nfreqs[1] = -1;
  Ndiffs = 2;					//number of diffs to be calculated (max 10)
  is[0] = offset+1;
  is[1] = -1;	

  continued_Nlambdas = 5;				//number of different lambdas to be calculated (max 10)
  continued_Nfreqs[0] = 1;				//the number of freqs used for calculation of eachj lambda
  continued_Nfreqs[1] = 2;
  continued_Nfreqs[2] = 5;
  continued_Nfreqs[3] = 20;
  continued_Nfreqs[4] = -1;
  continued_Ndiffs = 8;					//number of diffs to be calculated (max 10)
  continued_is[0] = 1;
  continued_is[1] = 2;			
  continued_is[2] = 3;			
  continued_is[3] = 4;			
  continued_is[4] = 5;			
  continued_is[5] = 6;			
  continued_is[6] = 7;						
  continued_is[7] = -1;	

  returned_Nfreq = 1;				// parameters for lambda that is the return value of CalculateLambda method. Number of freqs used and
  returned_FromContinued = true;				// whether input or continued X is used

  DoContinued = true;				// calculate and output  lambdas and diffs from continued X
  sprintf(OutputFN,"");				// output options
  DoOutput = true;				// if set to false, only returned lambda is calculated and no output files are created
  DoPrintOutContinuedX = false;
}

void LambdaCalculator::Initialize()
{ 
  for(int j=0; j<3; j++) 
  {
    rX[j] = new double[N];
    cXw[j] = new complex<double>[N];
    cXiw[j] = new complex<double>[M];
  }
  iw = new double[M];

  Initialized = true; 
}

void LambdaCalculator::ReleaseMemory()
{
  if (not Initialized) return;
  for(int j=0; j<3; j++) 
  {
    delete [] rX[j];
    delete [] cXw[j];
    delete [] cXiw[j];
  }
  delete [] iw;
}

//--------------------------------------------------------------//

void LambdaCalculator::SetN(int N)
{
  ReleaseMemory();
  this->N = N;
  Initialize();
  SetT(T);
}

void LambdaCalculator::SetT(double T)
{
  this->T = T;
  for(int m = 0; m < M; m++)
    iw[m] = (2.0*m+1.0)*pi*T;
}

void LambdaCalculator::SetNMT(int N, int M, double T)
{
  ReleaseMemory();
  this->N = N;
  this->M = M;
  Initialize();
  SetT(T);
  printf("LambdaCalculator: N: %d M: %d T: %.3f\n", this->N, this->M, this->T);
}

void LambdaCalculator::SetOmega(double* omega)
{
  this->omega = omega;
}

void LambdaCalculator::SetOffset(int offset)
{
  this->offset = offset;
  printf("LambdaCalculator: offset: %d\n", this->offset);
}

void LambdaCalculator::SetLambdasPrinted(int Nlambdas, const int* Nfreqs)
{
  this->Nlambdas = Nlambdas;
  for(int i=0; i < ((Nlambdas<=10) ? Nlambdas : 10) ; i++) this->Nfreqs[i] = Nfreqs[i];
}

void LambdaCalculator::SetContinuedLambdasPrinted(int continued_Nlambdas, const int* continued_Nfreqs)
{
  this->continued_Nlambdas = continued_Nlambdas;
  for(int i=0; i < ((continued_Nlambdas<=10) ? continued_Nlambdas : 10) ; i++) this->continued_Nfreqs[i] = continued_Nfreqs[i];
}

void LambdaCalculator::SetDiffsPrinted(int Ndiffs, const int* is)
{
  this->Ndiffs = Ndiffs;
  for(int i=0; i < ((Ndiffs<=10) ? Ndiffs : 10) ; i++) this->is[i] = is[i];
}

void LambdaCalculator::SetContinuedDiffsPrinted(int continued_Ndiffs, const int* continued_is)
{
  this->continued_Ndiffs = continued_Ndiffs;
  for(int i=0; i < ((continued_Ndiffs<=10) ? continued_Ndiffs : 10) ; i++) this->is[i] = is[i];
}


void LambdaCalculator::SetLambdaReturned(int returned_Nfreq, bool returned_FromContinued)
{
  this->returned_Nfreq = returned_Nfreq;
  this->returned_FromContinued = returned_FromContinued;
}

void LambdaCalculator::SetOutputFileName(const char* FN)
{
  sprintf(OutputFN,"%s",FN);
}

void LambdaCalculator::ResetCounter(int counter)
{
  this->counter = counter;
}

//--------------------------------------------------------//

void LambdaCalculator::ContinueToImagAxis(int N, complex<double>* Xw, double* w, int M, complex<double>* Xiw, double* iw, double T, const char* outputFN)
{
  complex<double>* x = new complex<double>[N];
  for(int m=0; m<M; m++)
  { 
    if (T>0) iw[m]=(2.0*m+1.0)*pi*T; // if T<0 use the values provided in iw
    for(int i=0; i<N; i++)
      x[i] = imag(Xw[i]) / complex<double>( w[i], iw[m] );
    
    Xiw[m] = -1.0/(pi)*TrapezIntegral(N, x, w);
  }
  delete [] x; 

  if (outputFN!=NULL)
    PrintFunc(outputFN, M, Xiw, iw); 
}


double LambdaCalculator::CalcLambda(int offset, int Nfreq, int N, complex<double>* X0, complex<double>* X1, complex<double>* X2)
{
  double num = 0.0;
  double denom = 0.0;
  for (int i=offset+1; i <= ((Nfreq>0) ? offset+Nfreq : N); i++)
  {
    num +=    abs(  ( X0[i] - X1[i] ) 
                   *( X1[i] - X2[i] ) );
    denom += sqr(abs( (X0[i] - X1[i]) ));
  }
  return 1.0 - num/denom;
}


double LambdaCalculator::CalcDiff(int offset, int i, int N, complex<double>* X1, complex<double>* X2 )
{
  if (i<0)
  { double max_diff=0.0;    
    for(int j = offset; j<N; j++)
    { double diff = abs(X2[j]-X1[j]);
      if (diff > max_diff) max_diff=diff;
    }
    return max_diff;
  }
  else
    return abs(X2[offset+i]-X1[offset+i]);
}

void  LambdaCalculator::PrintOut() 
{
  FILE* f; 
  char FN[300];
  if (counter>=2)
  { 
    sprintf(FN,"lambdas%s",OutputFN);
    if (counter==2) f = fopen(FN,"w"); else f = fopen(FN,"a");
    fprintf(f,"%d ",counter);
    for(int i = 0; i < Nlambdas; i++)
      fprintf(f,"%.15le ",lambdas[i]);
    fprintf(f,"\n"); 
    fclose(f);
  }

  if (counter>=1)
  {
    sprintf(FN,"diffs%s",OutputFN);
    if (counter==1) f = fopen(FN,"w"); else f = fopen(FN,"a");
    fprintf(f,"%d ",counter);
    for(int i = 0; i < Ndiffs; i++)
      fprintf(f,"%.15le ",diffs[i]);
    fprintf(f,"\n");
    fclose(f);
  }

  if (DoContinued)
  {
    if (counter>=2)
    { 
      sprintf(FN,"continued_lambdas%s",OutputFN);
      if (counter==2) f = fopen(FN,"w"); else f = fopen(FN,"a");
      fprintf(f,"%d ",counter);
      for(int i = 0; i < continued_Nlambdas; i++)
        fprintf(f,"%.15le ",continued_lambdas[i]);
      fprintf(f,"\n");
      fclose(f);
    }

    if (counter>=1)
    { 
      sprintf(FN,"continued_diffs%s",OutputFN);
      if (counter==1) f = fopen(FN,"w"); else f = fopen(FN,"a");
      fprintf(f,"%d ",counter);
      for(int i = 0; i < continued_Ndiffs; i++)
        fprintf(f,"%.15le ",continued_diffs[i]);
      fprintf(f,"\n");
      fclose(f);
    }

    if (DoPrintOutContinuedX)
    {
      sprintf(FN,"continued_X%s.it%d",OutputFN,counter);
      PrintFunc(FN,M,cXiw[2],iw);
    }
  }
}	

void  LambdaCalculator::ShiftAddAndContinue(complex<double>* cX)
{
  for(int i=0; i<N; i++)
  { cXw[0][i]=cXw[1][i];
    cXw[1][i]=cXw[2][i];
    cXw[2][i] = cX[i];
  } 
  if (DoContinued)
  { for(int i=0; i<M; i++)
    { cXiw[0][i]=cXiw[1][i];
      cXiw[1][i]=cXiw[2][i];
    }
    ContinueToImagAxis(N, cXw[2], omega, M, cXiw[2], iw);
  }
}

//===========================================================//

double LambdaCalculator::CalculateLambda(complex<double>* cX)
{
  ShiftAddAndContinue(cX);

  for(int i=0; i<Nlambdas; i++)
    lambdas[i] =  CalcLambda(offset, Nfreqs[i], N, cXw[0], cXw[1], cXw[2]);

  for(int i=0; i<Ndiffs; i++)
    diffs[i] =  CalcDiff( offset, is[i], N, cXw[1], cXw[2] );

  if (DoContinued)
  {
    for(int i=0; i<continued_Nlambdas; i++)
      continued_lambdas[i] =  CalcLambda(0, continued_Nfreqs[i], M, cXiw[0], cXiw[1], cXiw[2]);

    for(int i=0; i<continued_Ndiffs; i++)
      continued_diffs[i] =  CalcDiff(0, continued_is[i], M, cXiw[1], cXiw[2]);
  }

  PrintOut();

  counter++;

  return
    CalcLambda(  (returned_FromContinued) ? 0 : offset, 
                 returned_Nfreq, 
                 (returned_FromContinued) ? M : N, 
                 (returned_FromContinued) ? cXiw[0] : cXw[0], 
                 (returned_FromContinued) ? cXiw[1] : cXw[1], 
                 (returned_FromContinued) ? cXiw[2] : cXw[2] );
}

//===========================================================//



