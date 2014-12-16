#include "SIAM.h"
#include "routines.h"
#include "Broyden.h"
#include "GRID.h"
#include "Result.h"
#include "Input.h"

#ifdef _OMP
#include <omp.h>
#endif

//================== Constructors/DEstructors ====================//

void SIAM::Defaults()
{
  U = 2.0;
  T = 0.05;
  epsilon = 0;
  
  //broyden parameters
  MAX_ITS = 100; //default 100
  Accr = 1e-9; //default 1e-9

  // mu0 search
  max_tries = 2;
  UseBroydenFormu0 = true; 

  AmoebaScanStart = -2.0;
  AmoebaScanEnd = 2.0; 
  AmoebaScanStep = 0.2;
  AmoebaMaxIts = 100;
  AmoebaForceScanAndPrintOut = false;

  //broadening
  eta = 5e-2;
   
  //options
  CheckSpectralWeight = false; //default false
  UseMPT_Bs = false; //default false
  isBethe = false;
  mu0Fixed = false;
  mu0ismu = false;
  GfromDelta = false;

  UseLatticeSpecificG = false;
  t = 0.5;
  LatticeType = DOStypes::SemiCircle;
}

SIAM::SIAM()
{
  Defaults();
}

SIAM::SIAM(const char* ParamsFN)
{
  Defaults();
  this->ParamsFN.assign(ParamsFN);
  cout << "-- INFO -- SIAM: Params File Name set to:" << this->ParamsFN << endl;

  Input input(ParamsFN);

  input.ReadParam(U,"SIAM::U");
  input.ReadParam(T,"SIAM::T");
  input.ReadParam(epsilon,"SIAM::epsilon");
  input.ReadParam(eta,"SIAM::eta");
  input.ReadParam(MAX_ITS,"SIAM::MAX_ITS");
  input.ReadParam(Accr,"SIAM::Accr");
  input.ReadParam(CheckSpectralWeight, "SIAM::CheckSpectralWeight");
  input.ReadParam(UseMPT_Bs,"SIAM::UseMPT_Bs");
  input.ReadParam(isBethe,"SIAM::isBethe");
  input.ReadParam(mu0Fixed,"SIAM::mu0Fixed");
  input.ReadParam(mu0ismu,"SIAM::mu0ismu");
  Setmu0ismu(mu0ismu);
  input.ReadParam(GfromDelta,"SIAM::GfromDelta");
  input.ReadParam(UseLatticeSpecificG,"SIAM::UseLatticeSpecificG");
  input.ReadParam(UseBroydenFormu0,"SIAM::UseBroydenFormu0");
  input.ReadParam(max_tries,"SIAM::max_tries");
  input.ReadParam(AmoebaScanStart,"SIAM::AmoebaScanStart");
  input.ReadParam(AmoebaScanEnd,"SIAM::AmoebaScanEnd");
  input.ReadParam(AmoebaScanStep,"SIAM::AmoebaScanStep");
  input.ReadParam(AmoebaMaxIts,"SIAM::AmoebaMaxIts");
  input.ReadParam(AmoebaForceScanAndPrintOut,"SIAM::AmoebaForceScanAndPrintOut");
}

SIAM::~SIAM()
{

}

//========================= INITIALIZERS ===========================//

void SIAM::SetT(double T)
{
  this->T = T;
}

void SIAM::SetEpsilon(double epsilon)
{
  this->epsilon = epsilon;  
}


void SIAM::SetU(double U)
{
  this->U = U;
}

void SIAM::SetUTepsilon(double U, double T, double epsilon)
{
  SetU(U);
  SetT(T);
  SetEpsilon(epsilon);
}


void SIAM::SetBroydenParameters(int MAX_ITS, double Accr)
{
  this->MAX_ITS = MAX_ITS;
  this->Accr = Accr;
}

void SIAM::SetBroadening(double eta)
{
  this->eta = eta;
}

void SIAM::SetIsBethe(bool isBethe)
{
  this->isBethe = isBethe;
}

void SIAM::Setmu0Fixed(bool mu0Fixed)
{
  this->mu0Fixed = mu0Fixed;
}

void SIAM::Setmu0ismu(bool mu0ismu)
{
  this->mu0ismu = mu0ismu;
  if (mu0ismu) mu0Fixed=true;
}
void SIAM::SetUseLatticeSpecificG(bool UseLatticeSpecificG, double t, int LatticeType)
{
  this->UseLatticeSpecificG = UseLatticeSpecificG;
  this->t = t;
  this->LatticeType = LatticeType;
}

//========================= RUN SIAM EITH FIXED Mu ==========================//

bool SIAM::Run(Result* r) //output
{  
  this->r = r;
  N = r->grid->get_N();
  grid = r->grid;
  get_fermi();

  Clipped = false;
  
  if ((r->mu==0)&&(epsilon==-U/2.0)) 
    SymmetricCase = true;
  else 
    SymmetricCase = false;  
  
  printf("    -------%s SIAM: mu=%.3f, U=%.3f, T=%.3f, epsilon=%.3f -------\n", (SymmetricCase) ? "Symmetric" : "Asymmetric", r->mu, U, T, epsilon);
  
 
  //----- initial guess ------// 
  r->n = 0.5;  
  mu0 = r->mu0;  
  if ((!SymmetricCase)and(UseMPT_Bs))
     MPT_B = epsilon;
  else
     MPT_B = 0.0;

  complex<double>* V = new complex<double>[2];
  V[0] = mu0; 
  V[1] = MPT_B;
  //---------------------------//

  //------ SOLVE SIAM ---------//
  double mu0inits [] = {0.0, 1.0, -1.0, -0.8, 2.0, 1.5, -1.5,  2.5, -2.5, 
                        -2.0, 0.05, 0.8, 0.1, -0.1, 0.3, -0.3, 0.5, -0.5, -0.05, 
                        0.4, -0.4, 0.6, -0.6, 2.3, -2.3, 2.8, -2.8, 1.8, -1.8  }; 
  if ((SymmetricCase)or(mu0Fixed)) 
    //mu0 and n are known => there's no solving of system of equations
    SolveSiam(V);
  else 
  { bool failed = !UseBroydenFormu0;  
    int c = 0;
    while ( UseBroydenFormu0 and ( UseBroyden<SIAM>(2, MAX_ITS, Accr, &SIAM::SolveSiam, this, V) != 1 ) )
    { c++;
      if ( (c >= max_tries) or (c >= sizeof(mu0inits)/sizeof(double) - 1 ) )
      {
        printf("\n\n\n\n==== SIAM ERROR ====: Broyden mu0 search failed to converge. Now switching to amoeba ...\n\n\n\n");
        failed = true;  
        break;
      }
      V[0] = mu0inits[c]; 
      V[1] = MPT_B;
      printf("==================== ====================== ========== TRYING new mu0 int!!! c = %d, mu0init = %f\n\n\n",c, mu0inits[c]);
    };
    //use broyden to solve system of two equations
    if (failed) 
    {  V[0] = mu0inits[0]; 
       V[1] = MPT_B;
  
       Amoeba(Accr, V); 
    }
  }
  delete [] V;
  //----------------------------//

  //output spectral weight if opted
  if (CheckSpectralWeight)
  {
    printf("        Spectral weight G: %fe\n", -imag(TrapezIntegralMP(N, r->G, r->omega))/pi);
    printf("        Spectral weight G0: %fe\n", -imag(TrapezIntegralMP(N, r->G0, r->omega))/pi);
  }

  r->mu0 = mu0;
  printf("        SIAM: mu0 = %f, n = %f\n",r->mu0, r->n);


  #pragma omp parallel for
  for (int i=0; i<N; i++)
    r->DOS[i] = - imag(r->G[i]) / pi;

  return false;//Clipped;
}

//========================== RUN SIAM With FIXED n ==============================//
// applicable ONLY in solving Clean Hubbard Model which implies epsilon = 0 and NIDOS is needed on input.
//NOTE that MPT Bs will ALWAYS be one iteration late. They will converge to their real values
//when the DMFT loop converges. First DMFT Iteration is ALWAYS solved WITHOUT MPT Bs.

//TODO in case of asym NIDOS, mu and mu0 are not known EVEN FOR n=0.5 !!!!

bool SIAM::Run_CHM(Result* r) //output
{  
  this->r = r;
  N = r->grid->get_N();
  grid = r->grid;
  get_fermi();
  
  Clipped = false;
  
  epsilon = 0;

  if (r->n==0.5) HalfFilling = true; 
  else HalfFilling = false;
  
  printf("    ------- SIAM for CHM: n=%.3f, U=%.3f, T=%.3f, epsilon=%.3f -------\n", r->n, U, T, epsilon);
  
  if (HalfFilling) 
  {
    r->mu = 0.5*U;
    mu0 = 0.0;
    MPT_B = 0.0;
    MPT_B0 = 0.0;
    SymmetricCase = true;
  }
  else
  { mu0 = r->mu0;
    MPT_B = 0.0;
    MPT_B0 = 0.0;
    SymmetricCase = false;
  }

  complex<double>* V = new complex<double>[1];

  if (not mu0ismu)
  { //------initial guess---------//
    V[0] = mu0; //initial guess is always the last mu0. in first DMFT iteration it is 0
    //---------------------------//
    if (UseMPT_Bs)  printf("     MPT: B = %fe, B0 = %fe\n", MPT_B, MPT_B0);  

    //----------------- CALCULATION ----------------------//
    if ((HalfFilling)or(mu0Fixed))//and (SymmetricCase))
      get_G0();
    else
    { printf("         SIAM: about to calc mu0. at the moment: mu0 = %.3f mu=%.3f\n",r->mu0, r->mu);
      double initGuesses [] = {-1.0, 1.0, 0.3, -0.3, 0.1, -0.1, 
                               -0.8, 0.8, -0.6, 0.6, -0.7, 0.7,
                               -3.0, 3.0, 0.9, -0.9, 0.05, -0.05, 
                                0.5, -0.5, 0.2, -0.2, 2.0, -2.0};
    
      bool converged = false; 

      int i;
      for (i=0; ( (i<sizeof(initGuesses)/sizeof(double)) and (i<max_tries) ); i++)
      {  printf("------ SIAM: trying with init guess: %f\n",real(V[0]));
         converged = UseBroyden<SIAM>(1, 50, 1e-8, &SIAM::get_G0, this, V);  
         if (converged) break;
         else V[0] = initGuesses[i];
      }
      if ((i==max_tries)and(!converged))
      {  V[0] = initGuesses[0]; 
    
         Amoeba_CHM(Accr, V); 
      }
    }
    PrintFunc("G0", N, r->G0, r->omega);
    printf("    mu0 = %f\n", mu0);
  
    get_As();
    get_Ps();
    get_SOCSigma();
  }

  V[0] = r->mu;
  
  if (HalfFilling)//and (SymmetricCase))
  { if (isBethe)
    {  
      get_Sigma();
      get_G();
    }
    else
      get_G_CHM();
  }
  else
  { if (isBethe)
      UseBroyden<SIAM>(1, MAX_ITS, 1e-8, &SIAM::get_G, this, V);  
    else
    { 
      double initGuesses [] = {-1.0, 1.0, -0.8, 0.8, -0.6, 0.6, 0.2, -0.2, 2.0, -2.0};
    
      bool converged = false; 

      for (int i=0; i<sizeof(initGuesses)/sizeof(double); i++)
      {  printf("------ SIAM: trying with init guess: %f\n",real(V[0]));
         converged = UseBroyden<SIAM>(1, MAX_ITS, 1e-8, &SIAM::get_G_CHM, this, V);
         if ((converged)/*and(!Clipped)*/){ /*Clipped = false;*/ break; }
         else V[0] = initGuesses[i];
      }


    }

      
  }
  MPT_B = get_MPT_B();
  MPT_B0 = get_MPT_B0();

  printf("    mu = %f\n", r->mu);

  //delete [] V;
  //-----------------------------------------------------//

  //output spectral weight if optioned
  if (CheckSpectralWeight)
  {
    printf("    n0: %.6f\n", get_n(r->G0));
    printf("    n:  %.6f\n", get_n(r->G));
  }

  // fill in DOS
  #pragma omp parallel for
  for (int i=0; i<N; i++)
    r->DOS[i] = - imag(r->G[i]) / pi;

  r->mu0 = mu0;

  return false;
}

//=================================== FUNCTIONS ===================================//

double SIAM::get_fermi(int i)
{
  return 1.0 / ( 1.0 + exp( r->omega[i]/T ) );
}

void SIAM::get_fermi()
{  
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
  { //printf("tid: %d i: %d\n",omp_get_thread_num(),i);
    r->fermi[i] = get_fermi(i);
  }
}

void SIAM::get_G0()
{
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
    r->G0[i] = complex<double>(1.0)
               / ( complex<double>(r->omega[i] + mu0, eta)
                   - r->Delta[i] ); 

}

void SIAM::get_G0(complex<double>* V)
{
  mu0 = real(V[0]);

  get_G0();

  V[0] = mu0 + get_n(r->G0) - r->n;
  printf("get_G0: V[0] = %.5f\n",real(V[0]));
} 


double SIAM::get_n(complex<double> X[])
{
  double* g = new double[N];
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
    g[i]=-(1/pi)*imag(X[i])*r->fermi[i];
  
  double n = TrapezIntegralMP(N, g, r->omega);
  delete [] g;
  return n; 
}

void SIAM::get_As() 
{
  #pragma omp parallel for
  for (int i=0; i<N; i++)
  { //printf("tid: %d i: %d\n",omp_get_thread_num(),i);
    r->Ap[i] = -imag(r->G0[i]) * r->fermi[i] / pi;
    r->Am[i] = -imag(r->G0[i]) * (1.0 - r->fermi[i]) / pi;
  }
}

void SIAM::get_Ps()
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

      //PrintFunc("p1.w%.3f",N,p1,r->omega);
      //PrintFunc("p2",N,p1,r->omega);

      //get Ps by integrating                           
      r->P1[i] = pi * TrapezIntegral(N, p1[i], r->omega);
      r->P2[i] = pi * TrapezIntegral(N, p2[i], r->omega);

      delete [] p1[i];
      delete [] p2[i];
  }

  delete [] p1;
  delete [] p2;
}

/*
void SIAM::get_Ps()
{
  double** p1 = new double*[N];
  double** p2 = new double*[N];

  #pragma omp parallel for
  for (int i=N/2; i<N; i++) 
  {   //go over only positive frequences
      double* allomega = new double[2*N];
      p1[i] = new double[2*N];
      p2[i] = new double[2*N];

      int counter = 0;
      for (int j=0; j<N; j++)
        if(r->omega[j]<r->omega[i]/2.0)
        {
          allomega[counter] = r->omega[j];
          p1[i][counter]    = r->Am[j] * grid->interpl(r->Ap, r->omega[j] - r->omega[i]);
          p2[i][counter]    = r->Ap[j] * grid->interpl(r->Am, r->omega[j] - r->omega[i]);
          counter++;
        }
      for (int j=0; j<N; j++)
        if(r->omega[j]>=-r->omega[i]/2.0)
        {
          allomega[counter]=r->omega[j]+r->omega[i];
          p1[i][counter]    = r->Ap[j] * grid->interpl(r->Am, r->omega[j] + r->omega[i]);
          p2[i][counter]    = r->Am[j] * grid->interpl(r->Ap, r->omega[j] + r->omega[i]);
          counter++;
        }

      if (i%100==-1) //change this to zero to printout integrands
      { char p1FN[300];
        sprintf(p1FN,"p1.w%.3f",r->omega[i]); 
        PrintFunc(p1FN,counter,p1[i],allomega);
        char p2FN[300];
        sprintf(p2FN,"p2.w%.3f",r->omega[i]); 
        PrintFunc(p2FN,counter,p2[i],allomega);
      }

      //get Ps by integrating                           
      r->P1[i] = pi * TrapezIntegral(counter, p1[i], allomega);
      r->P2[i] = pi * TrapezIntegral(counter, p2[i], allomega);
      r->P1[N-1-i]=r->P2[i];
      r->P2[N-1-i]=r->P1[i];

      delete [] p1[i];
      delete [] p2[i];
      delete [] allomega;
  }

  delete [] p1;
  delete [] p2;
}
*/

void SIAM::get_SOCSigma()
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

double SIAM::get_MPT_B0()
{
  if (!UseMPT_Bs) return 0.0;
  
  complex<double>* b0 = new complex<double>[N]; //integrand function
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
    b0[i] = r->fermi[i] * r->Delta[i] * r->G0[i];
  
  double mpt_b0 = epsilon - 1.0  * (2.0 * r->n - 1.0) * imag(TrapezIntegralMP(N, b0, r->omega))
                           / ( pi * r->n * (1.0 - r->n) ) ;
  delete [] b0;
  return mpt_b0;
}

double SIAM::get_MPT_B()
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

double SIAM::get_b()
{ //we used mu0 as (mu0 - epsilon - U*n) in G0, now we're correcting that
  printf("         SIAM::get_b : MPT_B = %.3f, MPT_B0 = %.3f\n",MPT_B, MPT_B0);  
  if (!SymmetricCase)
    return ( (1.0 - 2.0 * r->n) * U - r->mu + (mu0 + epsilon + U * r->n) 
                             - MPT_B0 + MPT_B ) 
           /             ( r->n * (1.0 - r->n) * sqr(U) );
  else return 0;
}

void SIAM::get_Sigma()
{
 
  if (!SymmetricCase)
  { printf("going through asymmetric\n");
    double b = get_b();    
    #pragma omp parallel for
    for (int i=0; i<N; i++) 
      r->Sigma[i] =  U*r->n + r->SOCSigma[i] 
                              / ( 1.0 - b * r->SOCSigma[i] );
    
  }
  else
  { printf("going through symmetric\n");
    #pragma omp parallel for
    for (int i=0; i<N; i++) 
      r->Sigma[i] =  U * r->n + r->SOCSigma[i];
  }

}

//---------------- Get G -------------------------------//

void SIAM::get_G()
{
  if (UseLatticeSpecificG) 
    #pragma omp parallel for
    for (int i=0; i<N; i++) 
    { complex<double> com = r->omega[i] + r->mu - r->Sigma[i];
      r->G[i] = LS_get_G(LatticeType, t, com);
    }
  else
  {
  
  
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
  {    
    r->G[i] =  1.0
               / (r->omega[i] + r->mu - epsilon - r->Delta[i] - r->Sigma[i]) ;
    
    if (ClipOff(r->G[i])) Clipped = true;
  }
  
  if (Clipped) printf("    !!!!Clipping G!!!!\n");

  }
}

void SIAM::get_G(complex<double>* V)
{
  r->mu = real(V[0]);

  get_G();

  V[0] = r->mu + get_n(r->G) - r->n;
} 

//---------------- Get G for CHM -------------------------//

void SIAM::get_G_CHM()
{ 
  get_Sigma();   
  
  if (UseLatticeSpecificG) 
    #pragma omp parallel for
    for (int i=0; i<N; i++) 
    { complex<double> com = r->omega[i] + r->mu - r->Sigma[i];
      r->G[i] = LS_get_G(LatticeType, t, com);
    }
  else
  {

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

    if (ClipOff(r->G[i])) Clipped = true;
    delete [] g[i];
  }
  
  delete [] g;    
  
  if (Clipped) printf("    !!!!Clipping G!!!!\n");

  }
}

void SIAM::get_G_CHM(complex<double>* V)
{
  
  r->mu = real(V[0]);
  if(mu0ismu)
  { mu0 = r->mu;
    get_G0();
    get_As();
    get_Ps();
    get_SOCSigma();
  }

  get_G_CHM();

  V[0] = r->mu + get_n(r->G) - r->n;
  printf("during broyden mu=%.5f n(G)=%.5f\n",r->mu, get_n(r->G));
} 
//------------------------------------------------------//


void SIAM::SolveSiam(complex<double>* V)
{
  mu0 = real(V[0]);
  MPT_B = real(V[1]);

  //--------------------//
  get_G0();

  r->n0 = get_n(r->G0); 
  MPT_B0 = get_MPT_B0();  

  get_As();
  get_Ps();
  get_SOCSigma();
  
  r->n = r->n0; 

  if (GfromDelta)
  { get_Sigma();   
    get_G();
  }
  else
    get_G_CHM();
   
  r->n = get_n(r->G);
  //--------------------//

  V[0] = mu0 + (r->n - r->n0); //we need to satisfy (get_n(G) == n) and 
  V[1] = get_MPT_B();                //                (MPT_B == get_MPT_B())
}

void SIAM::Amoeba(double accr, complex<double>* V)
{
  //x here stands for mu0
  
  double x_start = AmoebaScanStart;
  double x_end   = AmoebaScanEnd;
  double x_step  = AmoebaScanStep;

  int sign_old=0;
  double x;
  bool found = false;
  int try_count = 0;
  double x_candidate;
  double x_best=0, diff_best=10e+100;
  
  while( (not found) and (try_count<1) ) 
  {
     FILE* ScanFile;  
     if (AmoebaForceScanAndPrintOut)
     {
       char ScanFN[50];
       sprintf(ScanFN, "scan.eps%.3f",epsilon);
       ScanFile = fopen(ScanFN,"w");
     }
     for(x=x_start; x<x_end; x+=x_step)
     {
       V[0] = x;
       SolveSiam(V);
       
       if (AmoebaForceScanAndPrintOut)
       {  char FN[50];
          sprintf(FN,"siam.eps%.3f.mu0_%.3f",epsilon, mu0);
          r->PrintResult(FN);
       }
      
       double x_res=real(V[0]);
       printf("         mu0: %.15f n(G0)-n(G): %.2le step: %.2le true diff: %.3f n(G): %.3f r->n: %.3f\n",x, x-x_res, x_step, get_n(r->G0)- get_n(r->G),get_n(r->G),r->n);

       if (AmoebaForceScanAndPrintOut)
         fprintf(ScanFile,"%.15le %.15le %.15le %.15le\n", x, x-x_res, r->n, r->n0);

       if (sign_old==0) 
       { sign_old = int_sign(x-x_res);
         continue;
       }

       int sign = int_sign(x - x_res);
       if (abs(x-x_res) < diff_best) { x_best=x; diff_best = abs(x-x_res); };
       if ((sign_old!=sign) and (not found))
       {  x_candidate = x-x_step;
          found = true; 
          if (not AmoebaForceScanAndPrintOut) break; 
       }
    }
    try_count++;
    if (not found) { x_start *=2.0; x_end *= 2.0; x_step *= 2.0; printf("              mu0 candidate NOT found! now scanning a wider range...\n"); }
    if (AmoebaForceScanAndPrintOut) fclose(ScanFile);
  } 
 
  
  if (not found)
  {  printf("              mu0 candidate NOT found! setting mu0 to to best choice: mu0_best: %f diff: %.2le\n",x_best,diff_best);
     V[0] = x_best;
     SolveSiam(V);
  }
  else
  {
    printf("              mu0 candidate found! proceeding with aomeba...\n");  
    x = x_candidate;
    x_step *= 0.5;
    x += x_step;

    bool converged = false;
    int it = 0;
    while( (not converged) and (it<=AmoebaMaxIts) )
    { it ++;
      V[0] = x;
      SolveSiam(V);
      double x_res=real(V[0]);
      printf("         it: %d mu0: %.15f n(G0)-n(G): %.2le step: %le n0: %.3f n: %.3f\n", it, x, x-x_res, x_step, get_n(r->G0), get_n(r->G));
      converged = ( abs(x-x_res) < accr );
      int sign = int_sign(x - x_res);
      if (sign_old==sign)
         x_step = abs(x_step);
      else
         x_step = -abs(x_step); 
      x_step *= 0.5;
      x += x_step;
    }
    if (converged) printf("          Amoeba: desired accuracy reached!\n");
  }
  printf("         --- Amoeba DONE ---\n");
}


void SIAM::Amoeba_CHM(double accr, complex<double>* V)
{
  //x here stands for mu0
  
  double x_start = AmoebaScanStart;
  double x_end   = AmoebaScanEnd;
  double x_step  = AmoebaScanStep;

  int sign_old=0;
  double x;
  bool found = false;
  int try_count = 0;
  double x_candidate;
  double x_best=0, diff_best=10e+100;
  
  while( (not found) and (try_count<1) ) 
  {
     FILE* ScanFile;  
     if (AmoebaForceScanAndPrintOut)
     {
       char ScanFN[50];
       sprintf(ScanFN, "scan.eps%.3f",epsilon);
       ScanFile = fopen(ScanFN,"w");
     }
     for(x=x_start; x<x_end; x+=x_step)
     {
       V[0] = x;
       get_G0(V);
       
       if (AmoebaForceScanAndPrintOut)
       {  char FN[50];
          sprintf(FN,"siam.eps%.3f.mu0_%.3f",epsilon, mu0);
          r->PrintResult(FN);
       }
      
       double x_res=real(V[0]);
       printf("         mu0: %.15f n(G0)-n(G): %.2le step: %.2le true diff: %.3f n(G0): %.3f r->n: %.3f\n",x, x-x_res, x_step, get_n(r->G0)- r->n, get_n(r->G0), r->n);

       if (AmoebaForceScanAndPrintOut)
         fprintf(ScanFile,"%.15le %.15le %.15le %.15le\n", x, x-x_res, r->n, r->n0);

       if (sign_old==0) 
       { sign_old = int_sign(x-x_res);
         continue;
       }

       int sign = int_sign(x - x_res);
       if (abs(x-x_res) < diff_best) { x_best=x; diff_best = abs(x-x_res); };
       if ((sign_old!=sign) and (not found))
       {  x_candidate = x-x_step;
          found = true; 
          if (not AmoebaForceScanAndPrintOut) break; 
       }
    }
    try_count++;
    if (not found) { x_start *=2.0; x_end *= 2.0; x_step *= 2.0; printf("              mu0 candidate NOT found! now scanning a wider range...\n"); }
    if (AmoebaForceScanAndPrintOut) fclose(ScanFile);
  } 
 
  
  if (not found)
  {  printf("              mu0 candidate NOT found! setting mu0 to to best choice: mu0_best: %f diff: %.2le\n",x_best,diff_best);
     V[0] = x_best;
     get_G0(V);
  }
  else
  {
    printf("              mu0 candidate found! proceeding with aomeba...\n");  
    x = x_candidate;
    x_step *= 0.5;
    x += x_step;

    bool converged = false;
    int it = 0;
    while( (not converged) and (it<=AmoebaMaxIts) )
    { it ++;
      V[0] = x;
      get_G0(V);
      double x_res=real(V[0]);
      printf("         it: %d mu0: %.15f n(G0)-n: %.2le step: %le n0: %.3f n: %.3f\n", it, x, x-x_res, x_step, get_n(r->G0), r->n);
      converged = ( abs(x-x_res) < accr );
      int sign = int_sign(x - x_res);
      if (sign_old==sign)
         x_step = abs(x_step);
      else
         x_step = -abs(x_step); 
      x_step *= 0.5;
      x += x_step;
    }
    if (converged) printf("          Amoeba: desired accuracy reached!\n");
  }
  printf("         --- Amoeba DONE ---\n");
}

//================================= ROUTINES =================================//

void SIAM::get_G_from_Sigma(Result* r)
{ 
  int N = r->grid->get_N();
  complex<double>** g = new complex<double>*[N];
  #pragma omp parallel for
  for (int i=0; i<N; i++) 
  {
      
    //treat integrand carefully 
    double D = 0.0;
    complex<double> LogTerm = 0.0;
    if (abs(imag(r->Sigma[i]))<0.1)
    {
      D = r->grid->interpl(r->NIDOS, r->mu + r->omega[i] - real(r->Sigma[i]));
      LogTerm = complex<double>(D, 0.0) * log( (r->mu + r->omega[i] - r->Sigma[i] + r->omega[N-1])
                                              /(r->mu + r->omega[i] - r->Sigma[i] - r->omega[N-1]) );
    }

    //create integrand array
    g[i] = new complex<double>[N];  
    for (int j=0; j<N; j++)
      g[i][j] = complex<double>(r->NIDOS[j] - D, 0.0) 
             / ( r->mu + r->omega[i] - r->omega[j] - r->Sigma[i] ); 

    /*if (i % 20 == 0) 
    { complex<double>* denom = new complex<double>[N];  
      for (int j=0; j<N; j++)
        denom[j] = 1.0 / ( r->mu + r->omega[i] - r->omega[j] - r->Sigma[i] ); 

      char denomFN[100];
      sprintf(denomFN,"denom.w%.3f",r->omega[i]);
      PrintFunc(denomFN, N, denom, r->omega);

      char integFN[100];
      sprintf(integFN,"integ.w%.3f",r->omega[i]);
      PrintFunc(integFN, N, g[i], r->omega);

      delete [] denom;
    }
  */
    //integrate to get G 
    r->G[i] = TrapezIntegral(N, g[i], r->omega) + LogTerm ; 

    delete [] g[i];
  }
  
  delete [] g;    


}



//-----------------------Miscellaneous---------------------------------//

bool SIAM::ClipOff(complex<double> &X)
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


double SIAM::MatsFreq(int m)
{
  return 2.0*pi*T*(m+0.5);
}

void SIAM::GetGfOnImagAxis(int M, complex<double> * G_out)
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
