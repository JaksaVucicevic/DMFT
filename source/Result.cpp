#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "Result.h"
#include "GRID.h"
#include "routines.h"
#include "IBZ.h"

Result::Result(GRID* grid)
{
  Initialize(grid);
}

Result::Result(const Result &result)
{
  Initialize(result.grid);
  CopyFrom(result);
}

Result::~Result()
{
  ReleaseMemory();
}

void Result::Reset()
{
  ReleaseMemory();
  Initialize(grid);
}

void Result::Reset(GRID* grid)
{
  ReleaseMemory();
  Initialize(grid);
}

void Result::Initialize(GRID* grid)
{
  this->grid = grid;
  
  int N = grid->get_N();
  
  omega = new double[N];
  grid->assign_omega(omega);

  Delta = new complex<double>[N];
  fermi = new double[N];
  G0 = new complex<double>[N];
  Ap = new double[N];
  Am = new double[N];
  P1 = new double[N];
  P2 = new double[N];
  SOCSigma = new complex<double>[N];
  Sigma = new complex<double>[N];
  G = new complex<double>[N];
  DOS = new double[N];
  NIDOS = new double[N];
  DOSmed = new double[N];

  n=0.0;
  n0=0.0;
  mu=0.0;
  mu0=0.0;
}

void Result::ReleaseMemory()
{
  delete [] omega;

  delete [] Delta;
  delete [] fermi;          
  delete [] G0;  
  delete [] Ap;          
  delete [] Am;
  delete [] P1;         
  delete [] P2;
  delete [] SOCSigma;
  delete [] Sigma;
  delete [] G;
  delete [] DOS;
  delete [] NIDOS;
  delete [] DOSmed;
}

void Result::PrintResult(const char* ResultFN)
{ 
  FILE *f;
  f = fopen(ResultFN, "w+");
  
  fprintf(f,"# n = %le n0 = %le mu = %le mu0=%le\n",n,n0,mu,mu0);   

  int N = grid->get_N();
  int i;
  for (i=0; i<N; i++)
  { 
     // loop through and store the numbers into the file
    fprintf(f, "%.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le %.15le\n", 
                   omega[i], fermi[i],					//1 2
                   real(Delta[i]), imag(Delta[i]),			//3 4
                   real(G0[i]), imag(G0[i]), 				//5 6
                   Ap[i], Am[i], P1[i], P2[i],				//7 8 9 10 
                   real(SOCSigma[i]), imag(SOCSigma[i]), 		//11 12
                   real(Sigma[i]), imag(Sigma[i]),			//13 14
                   real(G[i]), imag(G[i]),				//15 16
                   DOS[i], NIDOS[i], DOSmed[i]);			//17 18 19 
                   
                
  }
  fclose(f);
}

bool Result::ReadFromFile(const char* ResultFN)
{ 
  FILE *f;
  f = fopen(ResultFN, "r");
  if (f==NULL) { return false; }

  char rstLine[1000];
  fgets ( rstLine, 1000, f );

  char * pch;
  printf ("rstline: %s\n",rstLine);
  pch = strtok (rstLine,"=");
  int counter=1;
  while (pch != NULL)
  { 
    switch (counter)
    {
      case 2: n=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n); break;
      case 3: n0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,n0); break;
      case 4: mu=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu); break;
      case 5: mu0=atof(pch); printf ("%d: %s, atof: %le \n",counter,pch,mu0); break;
    }
          
    pch = strtok (NULL, "=");
    counter++; 
  }
  

  int N = grid->get_N();
  int i;
  for (i=0; i<N; i++)
  { double o, fer, rd, id, rg0, ig0, ap, am, p1, p2, rsocs, isocs, rs, is, rg, ig, dos, nidos, dosmed;
     // loop through and store the numbers into the file
    fscanf(f, "%le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n", 
                   &o, &fer,			//1 2
                   &rd, &id,			//3 4
                   &rg0, &ig0, 			//5 6
                   &ap, &am, &p1, &p2,		//7 8 9 10 
                   &rsocs, &isocs, 		//11 12
                   &rs, &is,			//13 14
                   &rg, &ig,			//15 16
                   &dos, &nidos, &dosmed);	//17 18 19 

    omega[i] = o;
    fermi[i] = fer;					
    Delta[i] = complex<double>(rd,id);			
    G0[i] = complex<double>(rg0,ig0); 				
    Ap[i] = ap; Am[i]=am; P1[i]=p1; P2[i]=p2;		
    SOCSigma[i] = complex<double>(rsocs,isocs); 	
    Sigma[i] = complex<double>(rs,is);			
    G[i] = complex<double>(rg,ig);			
    DOS[i] = dos; NIDOS[i] = nidos; DOSmed[i]=dosmed;	                                  
  }
  fclose(f);
  
  return true;
}


void Result::CopyFrom(const Result &result)
{
  Reset(result.grid);

  int N = result.grid->get_N();

  for (int i=0; i<N; i++)
  {
    //printf("m"); 
    omega[i] = result.omega[i];
    Delta[i] = result.Delta[i];
    fermi[i] = result.fermi[i];
    G0[i] = result.G0[i];
    Ap[i] = result.Ap[i];
    Am[i] = result.Am[i];
    P1[i] = result.P1[i];
    P2[i] = result.P2[i];
    SOCSigma[i] = result.SOCSigma[i];
    Sigma[i] = result.Sigma[i];
    G[i] = result.G[i];
    DOS[i] = result.DOS[i];
    NIDOS[i] = result.NIDOS[i];
    DOSmed[i] = result.DOSmed[i];
  }

  n = result.n;
  mu = result.mu;
  mu0 = result.mu0;
}

void Result::PrintModel(double U, double T)
{
  FILE *ModelFile, *InitImageFile;
  char FNModel[100], FNInitImage[100];
  sprintf(FNModel, "model.U%.3f.T%.3f", U, T);
  sprintf(FNInitImage, "initimage.U%.3f.T%.3f", U, T);
  ModelFile= fopen(FNModel,"w");
  InitImageFile= fopen(FNInitImage,"w");
          
  int nw = 200;  
  double dw = 0.025;
               
  for(int i=-nw; i<=nw; i++)
  {  fprintf(ModelFile, "%.15le  %.15le    %.15le\n", i*dw, -1.0/pi*imag(grid->interpl(G,i*dw)), dw ); 
     fprintf(InitImageFile, "%.15le   %.15le\n", i*dw, -1.0/pi*imag(grid->interpl(G,i*dw))); 
  }

  fclose(ModelFile);
  fclose(InitImageFile);
}

double Result::NIConductivity(double T, double mu, int Neps, int Nnu, const char * integrandFN)
{
  FILE* integrandFile;
  if (integrandFN!=NULL)
    integrandFile = fopen(integrandFN,"w");

  double sum = 0.0;
  double k = 12.0; 

  //double dnu = 2.0 * k * T / (double) Nnu;
  //for (double nu = -k*T; nu < k*T; nu += dnu )
 
  double nu_start = -3.0;
  double nu_end = 3.0;

  double dnu = (nu_end-nu_start)/ (double) Nnu;
 
  double eta=0.01; 

  for (double nu = nu_start; nu < nu_end; nu += dnu )
  { 
    double eps_start=-3.0; 
    double eps_end=3.0;

    double deps = (eps_end-eps_start) / (double) Neps ;
    for (double eps = eps_start; eps < eps_end; eps += deps )
    { 
      double v = (abs(eps)>=1.0) ? 0.0 : sqrt( (1.0-sqr(eps))/3.0 );
      double rho0 = grid->interpl(NIDOS,eps);
 
      complex<double> G = 1.0/( nu + mu - eps + ii*eta);
      double rho = - imag(G) / pi;
      double fprim = 1.0 / ( 4.0 * T * sqr( cosh( nu/(2.0*T) 
                                                ) 
                                          )   
                           );

      double integrand = rho0 * sqr(rho) * sqr(v) * fprim;
      sum += integrand * deps;
      if (integrandFN!=NULL) fprintf(integrandFile,"%.15le %.15le %.15le %.15le %.15le %.15le\n", eps, nu, integrand, rho, fprim, rho0);  
    }
    if (integrandFN!=NULL) fprintf(integrandFile,"\n");  
  }
  if (integrandFN!=NULL) fclose(integrandFile);

  return 2.0 * pi * sum * dnu ;
}



double Result::Conductivity(double T, double mu, int Neps, int Nnu, const char * integrandFN, bool excludeSmallOmega)
//--------------------- DC CONDUCTIVITY --------------------------// Neps,Nnu ~ 400 or 800
{
  FILE* integrandFile;
  if (integrandFN!=NULL)
    integrandFile = fopen(integrandFN,"w");

  double sum = 0.0;
  
  double k = 20.0;		//determines the nu range
  double W = 10.0;		//determines the epsilon range
  double NIDOS_EDGE = 1.0;	//edge of the non-interacting band

  //double dnu = 2.0 * k * T / (double) Nnu;
  //for (double nu = -k*T; nu < k*T; nu += dnu )
 
  double nu_start = -k*T;
  double nu_end = k*T;

  double dnu = (nu_end-nu_start)/ (double) Nnu;
 
  for (double nu = nu_start; nu < nu_end; nu += dnu )
  { 
    //---------- ONLY FOR INSULATOR -----------//
    if ((abs(nu)<0.1)and(excludeSmallOmega))
      continue;
    //-----------------------------------------//

    complex<double> Sigma_nu=grid->interpl(Sigma,nu);
    double eps_center = nu+mu-real(Sigma_nu);
    double eps_start=eps_center-W*abs(imag(Sigma_nu)); 
    double eps_end=eps_center+W*abs(imag(Sigma_nu));
    if (eps_start < -NIDOS_EDGE) eps_start = -1.0;
    if (eps_end > NIDOS_EDGE)    eps_end = 1.0;
    double deps = (eps_end-eps_start) / Neps ;
    for (double eps = eps_start; eps < eps_end; eps += deps )
    { 
      double v = (abs(eps)>=NIDOS_EDGE) ? 0.0 : 1.0;//sqrt( (1.0-sqr(eps))/3.0 );
      double rho0 = (abs(eps)>=NIDOS_EDGE) ? 0.0 : grid->interpl(NIDOS,eps);
 
      complex<double> G = 1.0/( nu + mu - eps - Sigma_nu);
      double rho = - imag(G) / pi;
      double fprim = 1.0 / ( 4.0 * T * sqr( cosh( nu/(2.0*T) 
                                                ) 
                                          )   
                           );

      double integrand = rho0 * sqr(rho) * sqr(v) * fprim;
      sum += integrand * deps;
      if (integrandFN!=NULL) fprintf(integrandFile,"%.15le %.15le %.15le %.15le %.15le %.15le\n", eps, nu, integrand, rho, fprim, rho0);  
    }
    if (integrandFN!=NULL) fprintf(integrandFile,"\n");  
  }
  if (integrandFN!=NULL) fclose(integrandFile);

  return 2.0 * pi * sum * dnu ;
}


double Result::TriangularConductivity(double T, int Nkx, int Nky, int Nnu, const char * integrandFN) //bool excludeSmallOmega)
//--------------------- DC CONDUCTIVITY --------------------------// Neps,Nnu ~ 400 or 800
{
  FILE* integrandFile;
  if (integrandFN!=NULL)
    integrandFile = fopen(integrandFN,"w");


  IBZ ibz(IBZtypes::TriangularLattice, Nkx, Nky );
  printf("-----ibz ready\n");
  complex<double> sum = 0.0;
  
  double k = 20.0;		//determines the nu range

  double nu_start = -k*T;
  double nu_end = k*T;

  double dnu = (nu_end-nu_start)/ (double) Nnu;
 
  for (double nu = nu_start; nu < nu_end; nu += dnu )
  { printf("-- nu: %.3f\n",nu);
    //---------- ONLY FOR INSULATOR -----------//
    //if ((abs(nu)<0.1)and(excludeSmallOmega))
    //  continue;
    //-----------------------------------------//

    complex<double> Sigma_nu=grid->interpl(Sigma,nu);
    double fprim = 1.0 / ( 4.0 * T * sqr( cosh( nu/(2.0*T) 
                                              ) 
                                        )   
                           );

    for(int i=0; i<Nkx; i++)
    for(int j=0; j<Nky; j++)
    { 
      double v = ibz.velocity[i][j];
 
      complex<double> G = 1.0/( nu + mu - ibz.epsilon[i][j] - Sigma_nu);
      double rho = - imag(G) / pi;

      ibz.summand[i][j] = sqr(rho) * sqr(v) * fprim;

      //if (integrandFN!=NULL) fprintf(integrandFile,"%.15le %.15le %.15le %.15le %.15le %.15le\n", eps, nu, integrand, rho, fprim, rho0);  
    }
    sum+=ibz.sum();

    if (integrandFN!=NULL) fprintf(integrandFile,"\n");  
  }
  if (integrandFN!=NULL) fclose(integrandFile);

  return 2.0 * pi * real(sum) * dnu ;
}

double Result::getIntegrand(double w, double T, double mu, double nu, double eps, complex<double> Sigma_nu,complex<double> Sigma_nu_plus_w)
{
      double v = 1.0;//sqrt( (1.0-sqr(eps))/3.0 );
      double rho0 = grid->interpl(NIDOS,eps);
 
      complex<double> G = 1.0/( nu + mu - eps - Sigma_nu );
      complex<double> Gw = 1.0/( nu + w + mu - eps - Sigma_nu_plus_w );
      double rho = - imag(G) / pi;
      double rhow = - imag(Gw) / pi;
      double fnu = 1.0 / (1.0 + exp(nu/T));
      double fnuw = 1.0 / (1.0 + exp((nu+w)/T));
      double integrand = rho0 * rho * rhow * sqr(v) * (fnu-fnuw)/w;
      return integrand;
}

double Result::Conductivity(double w, double T, double mu, double Neps, double Nnu, const char * integrandFN)
{
  double sum = 0.0;
  FILE* integrandFile;
  if (integrandFN!=NULL)
    integrandFile = fopen(integrandFN,"w");

  double k=20.0; 
  
  double nu_start = -w-k*T;
  double nu_end = k*T;
  double dnu = (nu_end-nu_start) / Nnu; 
  for (double nu = nu_start; nu < nu_end; nu += dnu )
  {  

    //---------- ONLY FOR INSULATOR -----------//
    //if((abs(nu)<0.1)or(abs(nu+w)<0.1))  
    //  continue;
    //-----------------------------------------//


    //-//
    complex<double> Sigma_nu=grid->interpl(Sigma,nu);
    complex<double> Sigma_nu_plus_w=grid->interpl(Sigma,nu+w);

    //-//
    double eps_center1 = nu+mu-real(Sigma_nu);
    double eps_start1 = eps_center1  - 65.0*abs(imag(Sigma_nu));
    double eps_end1   = eps_center1  + 65.0*abs(imag(Sigma_nu));

    if (eps_start1 < -1.0) eps_start1 = -1.0;
    if (eps_end1 > 1.0)    eps_end1 = 1.0;
    double deps1 = (eps_end1-eps_start1) / Neps ;

    //-//
    double eps_center2 = nu+w+mu-real(Sigma_nu_plus_w);
    double eps_start2 = eps_center2  - 65.0*abs(imag(Sigma_nu_plus_w));
    double eps_end2   = eps_center2  + 65.0*abs(imag(Sigma_nu_plus_w));

    if (eps_start2 < -1.0) eps_start2 = -1.0;
    if (eps_end2 > 1.0)    eps_end2 = 1.0;
    double deps2 = (eps_end2-eps_start2) / Neps ;

    //-// 
    double eps_start = (deps1<deps2) ? eps_start1 : eps_start2;
    double eps_end = (deps1<deps2) ? eps_end1 : eps_end2;            
    double deps = 2.0*( (deps1<deps2) ? deps1 : deps2 );
    for (double eps = eps_start; eps < eps_end; eps += deps )
    { 
      double integrand = getIntegrand( w, T, mu, nu, eps, Sigma_nu, Sigma_nu_plus_w );
      sum += integrand * deps;
      if (integrandFN!=NULL) fprintf(integrandFile,"%.15le %.15le %.15le\n", eps, nu, integrand /*rho, (fnu-fnuw)/w*/);  
    }

    //-// 
    eps_start = (deps1>=deps2) ? eps_start1 : eps_start2;
    eps_end = (deps1>=deps2) ? eps_end1 : eps_end2;            
    deps = 2.0*( (deps1>=deps2) ? deps1 : deps2 );
    for (double eps = eps_start; eps < eps_end; eps += deps )
    { 
      if ( ( eps > ((deps1<deps2) ? eps_start1 : eps_start2)
           ) 
            and 
           ( eps < ((deps1<deps2) ? eps_end1 : eps_end2)
           ) 
         )
        continue;
      double integrand = getIntegrand( w, T, mu, nu, eps, Sigma_nu, Sigma_nu_plus_w );
      sum += integrand * deps;
      if (integrandFN!=NULL) fprintf(integrandFile,"%.15le %.15le %.15le\n", eps, nu, integrand /*rho, (fnu-fnuw)/w*/);  
    }

    //-//
    if (integrandFN!=NULL) fprintf(integrandFile,"\n");  
  }
  if (integrandFN!=NULL) fclose(integrandFile);

  return 2.0 * pi * sum * dnu ;
}


/*
double Result::Conductivity(double w, double T, double mu, double Neps, double Nnu, double Nscan)
{
//  double n = 100.0;
//  double k = 10.0;
  double sum = 0.0;
  //char integrandFN[100];
  //sprintf(integrandFN, "integrand.U%.3f.T%.3f",2.0*mu,T);
  //FILE* integrandFile = fopen(integrandFN,"w");
  double k=5.0;
  double nu = k*T/2.0;
  double eps_max = 0.0;
  double rho_max = 0.0;
  for (double eps = -1.0; eps < 1.0; eps += 2.0 / Nscan )
  {
    complex<double> G = 1.0/( nu + mu - eps - grid->interpl(Sigma,nu) );
    double rho = - imag(G) / pi;
    
    if (rho_max < rho) { rho_max = rho; eps_max = eps; }
  }
  double halfWidth=0;
  for (double eps = -1.0; eps < 1.0; eps += 2.0 / Nscan )
  {
    complex<double> G = 1.0/( nu + mu - eps - grid->interpl(Sigma,nu) );
    double rho = - imag(G) / pi;
    
    if (rho_max < 100*rho) { halfWidth = abs(eps_max-eps); break; } 
  }

  double angle = nu/eps_max;

  if (angle<0) k = 12.0; 
  printf("rho_max: %f eps_max: %f angle: %f halfWidth: %f\n",rho_max, eps_max, angle, halfWidth);

  k=10.0; 
  double dnu = 2.0 * k * T / Nnu;
  for (double nu = -k*T; nu < k*T; nu += dnu )
  { 
    double eps_start = -1.0;
    double eps_end = 1.0;
    //double eps_start = nu/angle - halfWidth;
    //double eps_end = nu/angle + halfWidth;
    //if ((eps_start < -1.0)or(angle<0)) eps_start = -1.0;
    //if ((eps_end > 1.0)or(angle<0)) eps_end = 1.0;
    double deps = (eps_end-eps_start) / Neps ;
    for (double eps = eps_start; eps < eps_end; eps += deps )
    { 
      double v = sqrt( (1.0-sqr(eps))/3.0 );
      double rho0 = grid->interpl(NIDOS,eps);
 
      complex<double> G = 1.0/( nu + mu - eps - grid->interpl(Sigma,nu) );
      complex<double> Gw = 1.0/( nu + w + mu - eps - grid->interpl(Sigma,nu) );
      double rho = - imag(G) / pi;
      double rhow = - imag(Gw) / pi;
      double fnu = 1.0 / (1.0 + exp(nu/T));
      double fnuw = 1.0 / (1.0 + exp((nu+w)/T));
      double integrand = rho0 * rho * rhow * sqr(v) * (fnu-fnuw)/w;
      sum += integrand * deps;
      //fprintf(integrandFile,"%.15le %.15le %.15le %.15le %.15le\n", eps, nu, integrand, rho, (fnu-fnuw)/w);  
    }
    //fprintf(integrandFile,"\n");  
  }
  //fclose(integrandFile);

  return 2.0 * pi * sum * dnu ;
}*/




void Result::ChargeSusceptibility(double T, double &chi1, double &chi3)
{
  int N = grid->get_N();
  double* intg = new double[N];
  for (int i=0; i<N; i++)
  {
    intg[i] = -(1.0/pi) * fermi[i] * imag(G[i]);
  }
  double n1 = TrapezIntegral(N,intg,omega);
 

  double dmu = 0.01;
  for (int i=0; i<N; i++)
  {
    double f = 1.0/(1.0 + exp( (omega[i] - dmu)/T) );
    intg[i] = -(1.0/pi) * f * imag(G[i]);
  }
  double n2 = TrapezIntegral(N,intg,omega);
  chi1 = (n2-n1)/dmu;
/*
  double* g = new double[N];
  for (int i=0; i<N; i++) // omega integral
  { 
    for (int j=0; j<N; j++) //epsilon integral
    {
      intg[j] =   imag(Sigma[i]) 
                * ( omega[i] + mu - omega[j] - real(Sigma[i]) ) 
                * NIDOS[j]
                / sqr(     sqr( omega[i] + mu - omega[j] - real(Sigma[i])  ) 
                        +  sqr(imag(Sigma[i])) 
                     ) 
                 ;
    }
    char fn[100];
    sprintf(fn,"intg.w%.3f",omega[i]);
    if (i%100==0) PrintFunc(fn,N,intg,omega);
    g[i] = fermi[i] * TrapezIntegral(N,intg,omega);   
  }
  PrintFunc("g",N,g,omega);
  chi3 = (2.0/pi) * TrapezIntegral(N,g,omega);

  delete [] g;
  printf("n1: %f n2: %f chi1: %f  chi3: %f \n",n1,n2,chi1,chi3);
*/

  dmu = 0.02;
  for (int i=0; i<N; i++) 
  {
      
    //treat integrand carefully 
    double D = 0.0;
    complex<double> LogTerm = 0.0;
    if (abs(imag(Sigma[i]))<0.1) 
    {
      D = grid->interpl(NIDOS, dmu + omega[i] - real(SOCSigma[i]));
      LogTerm = complex<double>(D, 0.0) * log( (dmu + omega[i] - SOCSigma[i] + omega[N-1])
                                              /(dmu + omega[i] - SOCSigma[i] - omega[N-1]) );
    }

    //create integrand array
    complex<double>* g = new complex<double>[N];  
    for (int j=0; j<N; j++)
      g[j] = complex<double>(NIDOS[j] - D, 0.0) 
             / ( dmu + omega[i] - omega[j] - SOCSigma[i] ); 
    
  
    //integrate to get G 
    G[i] = TrapezIntegral(N, g, omega) + LogTerm ; 

    delete [] g;
  }
  char fn[100];
  sprintf(fn,"NEW3/Gprobno.U%.3f.T%.3f", 2.0 * (real(Sigma[0])-real(SOCSigma[0])), T );
  PrintFunc(fn,N,G,omega);
  for (int i=0; i<N; i++)
    intg[i] = -(1.0/pi) * fermi[i] * imag(G[i]);
//      intg[i] = -(1.0/pi) * imag(G[i])
 //                         * (1.0/(1.0 + exp( (omega[i] - dmu)/T 
  //                                         ) 
    //                              )
      //                      ) ;
  double n3 = TrapezIntegral(N,intg,omega);
  chi3 = (n3-n1)/dmu;

  printf("n1: %f n2: %f n3: %f chi1: %f  chi3: %f mu: %f \n",n1,n2,n3,chi1,chi3,mu);
  delete [] intg;

}

void Result::PrintSpectralFunction(const char* FN, double (*eps)(double, double))
{
  double Kx = 4.0*pi/3.0;
  double ks [3][2] = {{0,0},{Kx,0},{Kx*3.0/4.0, Kx*sqrt(3.0)/4.0}};

  int N = 100;
  int Ns [3] = {N, N/2, (int)(N*sqrt(3.0)/2.0)};

  FILE* AFile = fopen(FN,"w");
  for(double w=-3.0; w<3.0; w+=0.01)
  { complex<double> Sigma_w=grid->interpl(Sigma,w);
    int counter = 0;
    for(int l=0; l<3; l++)
    { double dkx = ( ks[(l+1==3)?0:(l+1)][0] - ks[l][0] )/((double) Ns[l]);
      double dky = ( ks[(l+1==3)?0:(l+1)][1] - ks[l][1] )/((double) Ns[l]);
      for(int i=0; i<Ns[l]; i++)
      { 
        double kx = ks[l][0]+i*dkx;
        double ky = ks[l][1]+i*dky;
  
        double A = -(1.0/pi)*imag(1.0/( w + mu - eps(kx,ky) - Sigma_w ));
        fprintf(AFile,"%d %.15le %.15le\n", counter, w, A);
        counter++;
      }
    }  
    fprintf(AFile,"\n");  
  }
  fclose(AFile);
}









void Result::PrintOnImagAxis(complex<double> * X, int M, double T, const char* FN)
{ 
  int N = grid->get_N();
  complex<double>* g = new complex<double>[N];
  complex<double>* IAX = new complex<double>[M];
  double* mf = new double[M];

  for(int m=0; m<M; m++)
  { 
    mf[m]= 2.0*pi*T*(m+0.5);
    for(int i=0; i<N; i++)
      g[i] = imag(X[i]) / complex<double>( -omega[i], mf[m] );
    
    IAX[m] = -1/(pi)*TrapezIntegral(N, g, omega);
  }

  PrintFunc(FN, M, IAX, mf);

  delete [] IAX;
  delete [] g;
  delete [] mf;
}

void Result::PrintCondOnImagAxis(int N, double* X, double* w, int M, double T, const char* FN)
{ 
  double* g = new double[N];
  double* IAX = new double[M];
  double* mf = new double[M];

  for(int m=0; m<M; m++)
  { 
    mf[m]= 2.0*pi*T*(m+1);
    for(int i=0; i<N; i++)
      g[i] =            X[i] * mf[m] 
              /  (  sqr(w[i]) + sqr(mf[m])  );
    
    IAX[m] = (2.0/pi)*TrapezIntegral(N, g, w);
  }

  PrintFunc(FN, M, IAX, mf);

  delete [] IAX;
  delete [] g;
  delete [] mf;
}

double Result::get_ImGiw1(double T)
{ 
  int N = grid->get_N();
  complex<double>* g = new complex<double>[N];
  
  double w1 = pi*T;
  for(int i=0; i<N; i++)
    g[i] = imag(G[i]) / complex<double>( -omega[i], w1 );
    
  complex<double> Giw1 = -1.0/(pi)*TrapezIntegral(N, g, omega);  
  delete [] g;

  return imag(Giw1);
}



