#include "SIAM2.h"
#include "routines.h"
#include "Broyden.h"
#include "GRID.h"
#include "Result.h"
#include "Input.h"

#ifdef _OMP
#include <omp.h>
#endif

//================== Constructors/DEstructors ====================//

void SIAM2::Defaults()
{
  U = 2.0;
  T = 0.05;
  epsilon = 0;
  
  //broyden parameters
  MAX_ITS = 100; //default 100
  Accr = 1e-9; //default 1e-9

  //broadening
  eta = 5e-2;
   
  //options
  CheckSpectralWeight = false; //default false
  UseMPT_Bs = false; //default false

  t = 0.5;
  LatticeType = DOStypes::SemiCircle;

  FORMULA_FOR_G = FormulasForG::Integral;
  HalfFilling = false;
}

SIAM2::SIAM2()
{
  Defaults();
}

SIAM2::SIAM2(const char* ParamsFN)
{
  Defaults();
  this->ParamsFN.assign(ParamsFN);
  cout << "-- INFO -- SIAM2: Params File Name set to:" << this->ParamsFN << endl;

  Input input(ParamsFN);

  input.ReadParam(U,"SIAM2::U");
  input.ReadParam(T,"SIAM2::T");
  input.ReadParam(epsilon,"SIAM2::epsilon");
  input.ReadParam(eta,"SIAM2::eta");
  input.ReadParam(MAX_ITS,"SIAM2::MAX_ITS");
  input.ReadParam(Accr,"SIAM2::Accr");
  input.ReadParam(CheckSpectralWeight, "SIAM2::CheckSpectralWeight");
  input.ReadParam(UseMPT_Bs,"SIAM2::UseMPT_Bs");
}

SIAM2::~SIAM2()
{

}

//========================= INITIALIZERS ===========================//

void SIAM2::SetT(double T)
{
  this->T = T;
}

void SIAM2::SetEpsilon(double epsilon)
{
  this->epsilon = epsilon;  
}


void SIAM2::SetU(double U)
{
  this->U = U;
}

void SIAM2::SetUTepsilon(double U, double T, double epsilon)
{
  SetU(U);
  SetT(T);
  SetEpsilon(epsilon);
}


void SIAM2::SetBroydenParameters(int MAX_ITS, double Accr)
{
  this->MAX_ITS = MAX_ITS;
  this->Accr = Accr;
}

void SIAM2::SetBroadening(double eta)
{
  this->eta = eta;
}


void SIAM2::SetUseLatticeSpecificG(double t, int LatticeType)
{
  FORMULA_FOR_G = FormulasForG::LatticeSpecific;
  this->t = t;
  this->LatticeType = LatticeType;
}

//========================= RUN SIAM2 EITH FIXED Mu ==========================//

bool SIAM2::Run(Result* r) //output
{  
  this->r = r;
  N = r->grid->get_N();
  grid = r->grid;
  get_fermi();

  Clipped = false;
  
  printf("    ------- SIAM2: U=%.3f, T=%.3f, epsilon=%.3f %s-------\n", U, T, epsilon, (HalfFilling) ? "at Half filling" : "");
  printf("    -- params from result: n=%.3f, n0=%.3f, mu=%.3f mu0=%.3f -------\n", r->n, r->n0, r->mu, r->mu0);
 


  //------ SOLVE SIAM2 ---------//
  get_G0();

  get_As();
  get_Ps();
  get_SOCSigma();

  if (!HalfFilling) r->n0 = get_n(r->G0);
  printf(" n0 new: %f\n",r->n0);
  MPT_B0 = get_MPT_B0();  

  get_Sigma();   
  get_G();

  if (!HalfFilling) r->n = get_n(r->G);
  printf(" n new: %f\n",r->n);
  //----------------------------//
  
  //output spectral weight if optioned
  if (CheckSpectralWeight)
  {
    printf("        Spectral weight G: %fe\n",-(1.0/pi)*imag(TrapezIntegral(N,r->G,r->omega)));
    printf("        Spectral weight G0: %fe\n",-(1.0/pi)*imag(TrapezIntegral(N,r->G0,r->omega)));
  }

  //return Clipped;
  return false;
}

//=================================== FUNCTIONS ===================================//

double SIAM2::get_fermi(int i)
{
  return 1.0 / ( 1.0 + exp( r->omega[i]/T ) );
}

void SIAM2::get_fermi()
{  
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
  { //printf("tid: %d i: %d\n",omp_get_thread_num(),i);
    r->fermi[i] = get_fermi(i);
  }
}

void SIAM2::get_G0()
{
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
    r->G0[i] = complex<double>(1.0)
               / ( complex<double>(r->omega[i] + r->mu0 - epsilon - U*r->n, eta)
                   - r->Delta[i] ); 

}


double SIAM2::get_n(complex<double> X[])
{
  double* g = new double[N];
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
    g[i]=-(1/pi)*imag(X[i])*r->fermi[i];
  
  double n = TrapezIntegralMP(N, g, r->omega);
  delete [] g;
  return n; 
}

void SIAM2::get_As() 
{
  #pragma omp parallel for
  for (int i=0; i<N; i++)
  { //printf("tid: %d i: %d\n",omp_get_thread_num(),i);
    r->Ap[i] = -imag(r->G0[i]) * r->fermi[i] / pi;
    r->Am[i] = -imag(r->G0[i]) * (1.0 - r->fermi[i]) / pi;
  }
}

void SIAM2::get_Ps()
{
  double** p1 = new double*[N];
  double** p2 = new double*[N];

  #pragma omp parallel for
  for (int i=0; i<N; i++) 
  { 
      p1[i] = new double[N];
      p2[i] = new double[N];
      for (int j=0; j<N; j++)
      {  
         p1[i][j] = r->Am[j] * grid->interpl(r->Ap, r->omega[j] - r->omega[i]);
         p2[i][j] = r->Ap[j] * grid->interpl(r->Am, r->omega[j] - r->omega[i]);
      }

      //get Ps by integrating                           
      r->P1[i] = pi * TrapezIntegral(N, p1[i], r->omega);
      r->P2[i] = pi * TrapezIntegral(N, p2[i], r->omega);

      delete [] p1[i];
      delete [] p2[i];
  }

  delete [] p1;
  delete [] p2;
}

void SIAM2::get_SOCSigma()
{
    double** s = new double*[N];
    #pragma omp parallel for 
    for (int i=0; i<N; i++) 
    { //printf("tid: %d i: %d\n",omp_get_thread_num(),i);
      s[i] = new double[N];
      for (int j=0; j<N; j++) 
      {  //printf("tid: %d j: %d\n",omp_get_thread_num());
         s[i][j] =   grid->interpl(r->Ap, r->omega[i] - r->omega[j]) * r->P2[j] 
                   + grid->interpl(r->Am, r->omega[i] - r->omega[j]) * r->P1[j];
      }
                         
      //integrate
  
      r->SOCSigma[i] = complex<double>(0.0, - U*U * TrapezIntegral(N, s[i], r->omega) );    
     
      if (ClipOff( r->SOCSigma[i] )) Clipped = true ;
      delete [] s[i];
    }
    delete [] s;
  //int i;
  //cin >> i;
  if (Clipped) printf("    !!!Clipping SOCSigma!!!!\n");
  grid->KramarsKronig( r->SOCSigma );
}

double SIAM2::get_MPT_B0()
{
  if (!UseMPT_Bs) return 0.0;
  
  complex<double>* b0 = new complex<double>[N]; //integrand function
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
    b0[i] = r->fermi[i] * r->Delta[i] * r->G0[i];
  
  double mpt_b0 = epsilon - 1.0  * (2.0 * r->n0 - 1.0) * imag(TrapezIntegralMP(N, b0, r->omega))
                           / ( pi * r->n0 * (1.0 - r->n0) ) ;
  delete [] b0;
  return mpt_b0;
}

double SIAM2::get_MPT_B()
{
  if (!UseMPT_Bs) return 0.0;
  
  complex<double>* b = new complex<double>[N];
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
    b[i] = r->fermi[i] * r->Delta[i] * r->G[i]
           * ( (2.0 / U) * r->Sigma[i] - 1.0 );
  
  double mpt_b = epsilon - 1.0/( pi * r->n * (1.0 - r->n) ) 
                           * imag(TrapezIntegralMP(N, b, r->omega));
  delete [] b;
  return mpt_b;
}

double SIAM2::get_a()
{ 
  if (HalfFilling) return 1.0;
  return   (  r->n * (1.0 - r->n)  )
         / ( r->n0 * (1.0 - r->n0) );
}

double SIAM2::get_b()
{ 
  if (HalfFilling) return 0.0;
  return ( (1.0 - 2.0 * r->n) * U - MPT_B0 + MPT_B ) 
         /   ( r->n0 * (1.0 - r->n0) * sqr(U) );
}

void SIAM2::get_Sigma()
{
   double b = get_b();    
   double a = get_a();    
   printf("---get_Sigma: a: %.3f b: %.3f n: %.3f n0: %.3f\n",a,b,r->n,r->n0);
   #pragma omp parallel for
   for (int i=0; i<N; i++) 
     r->Sigma[i] =  U*r->n +       a * r->SOCSigma[i] 
                             / ( 1.0 - b * r->SOCSigma[i] );    
}

//---------------- Get G -------------------------------//

void SIAM2::get_G()
{
  switch(FORMULA_FOR_G) 
  { case 0:  
      printf("Lattice specific\n");
      #pragma omp parallel for
      for (int i=0; i<N; i++) 
      { complex<double> com = r->omega[i] + r->mu - r->Sigma[i];
        r->G[i] = LS_get_G(LatticeType, t, com);
      }
      break;
    case 1:   
      printf("No lattice\n");
      #pragma omp parallel for
      for (int i=0; i<N; i++) 
      {    
        r->G[i] =  1.0
                   / (r->omega[i] + r->mu - epsilon - r->Delta[i] - r->Sigma[i]) ;
      }
      break;
    case 2:
      printf("Integral\n");
      complex<double>** g = new complex<double>*[N];
      #pragma omp parallel for
      for (int i=0; i<N; i++) 
      {
      
        //treat integrand carefully 
        double D = 0.0;
        complex<double> LogTerm = 0.0;
        if (abs(imag(r->Sigma[i]))<0.1) 
        {
          D = grid->interpl(r->NIDOS, r->mu + r->omega[i] - real(r->Sigma[i]));
          LogTerm = complex<double>(D, 0.0) * log( (r->mu + r->omega[i] - r->Sigma[i] + r->omega[N-1])
                                                  /(r->mu + r->omega[i] - r->Sigma[i] - r->omega[N-1]) );
        }

        //create integrand array
        g[i] = new complex<double>[N];  
        for (int j=0; j<N; j++)
          g[i][j] = complex<double>(r->NIDOS[j] - D, 0.0) 
                    / ( r->mu + r->omega[i] - r->omega[j] - r->Sigma[i] ); 
    
  
        //integrate to get G 
        r->G[i] = TrapezIntegral(N, g[i], r->omega) + LogTerm ; 
 
        delete [] g[i];
      }
      delete [] g;    
      break;
    default: printf("Not implemented!!!\n"); exit(1);   
  }
}



//================================= ROUTINES =================================//

//-----------------------Miscellaneous---------------------------------//

bool SIAM2::ClipOff(complex<double> &X)
{
  if (imag(X)>0) 
  {
    X = complex<double>(real(X),-1e-5);
    return true;
  }
  else
    return false;
}

//------------------------ IMAG Axis ---------------------------//


double SIAM2::MatsFreq(int m)
{
  return 2.0*pi*T*(m+0.5);
}

void SIAM2::GetGfOnImagAxis(int M, complex<double> * G_out)
{ 
  complex<double>* g = new complex<double>[N];
  double* mf = new double[M];
  for(int m=0; m<M; m++)
  { 
    mf[m]=MatsFreq(m);
    for(int i=0; i<N; i++)
      g[i] = imag(r->G[i]) / complex<double>( -r->omega[i], mf[m] );
    
    G_out[m] = -1/(pi)*TrapezIntegral(N, g, r->omega);
  }
  delete [] g;
  delete [] mf;
}
