#include "MHCHM.h"
#include <complex>
#include "routines.h"
#include "GRID.h"
#include "Result.h"

using namespace std;

void MHCHM::Defaults()
{
  UseBethe = true;  
}

MHCHM::MHCHM() : CHM()
{
  Defaults();
}

MHCHM::MHCHM(const char* ParamsFN) : CHM(ParamsFN)
{
  Defaults();
  
}

void MHCHM::ReleaseMemory()
{
  printf("MHCHM release\n");
}

MHCHM::~MHCHM() 
{
  ReleaseMemory();
}

void MHCHM::SetUseBethe(bool UseBethe)
{
  
}

bool MHCHM::SolveSIAM()
{
  printf("entered...\n");
 
  rho = new complex<double>[N];
  rhot = new complex<double>[N];
  Sigmat = new complex<double>[N]; 

  get_G0();
  printf("got G0...\n");
  
  get_rhot();
  printf("got rhot...\n");

  get_Sigmat();
  printf("got Sigmat...\n");

  get_Sigma();
  printf("got Sigma...\n");

  get_G();
  printf("got G...\n");

  delete [] rho;
  delete [] rhot;
  delete [] Sigmat;

  return false;
}

void MHCHM::get_G0()
{
  for (int i = 0; i < N; i++)
    r->G0[i] = 1.0 
               / ( r->omega[i] - r->Delta[i] + ii*SIAMeta ) ;
}

void MHCHM::get_rhot()
{
  double domega = grid->get_domega();

  for(int i=0; i<N; i++) 
    rho[i] = - (1.0 / pi) * imag(r->G0[i]); 
 
  FastFourier(rho, rhot, domega, -1);
      
  PrintFunc("rhot", N, rhot);

}

void MHCHM::get_Sigmat()
{
  for (int j = 0; j<N/2; j++)
    Sigmat[j] = pow(U, 2.0) * pow(rhot[j], 3.0);

  // erase the upper half of times
  for (int j = N/2; j<N; j++)
    Sigmat[j] = 0.0;

  PrintFunc("Sigmat", N, Sigmat);
}

void MHCHM::get_Sigma()
{
  double dt = 2.0 * pi / (N*grid->get_domega());

  FastFourier(Sigmat, r->Sigma, dt, 1);
  
  bool Clipped;
  for(int i=0; i<N/2; i++)
  {   
    r->Sigma[i] *= - 2.0 * ii;
    
    // Clip Off Sigma
    if (imag(r->Sigma[i]) > 0.0) 
    {  r->Sigma[i] = complex<double>(real(r->Sigma[i]),0.0);
       Clipped = true;
    }
    
  }
  if (Clipped) printf("--- WARNING --- Sigma Clipped !!!\n");
  
  // erase the upper half of frequences
  for(int i=N/2; i<N; i++)
    r->Sigma[i] = 0.0;
}

void MHCHM::get_G()
{
  complex<double>* comega = new complex<double>[N];
  for (int i = 0; i<N; i++)
   comega[i] = r->omega[i] - r->Sigma[i] + ii*SIAMeta;

  InitG(DOStypes::SemiCircle, t, N, comega, r->G);  
}


//----------Fast Fourier-----------------------------//
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void MHCHM::FastFourier(complex<double>* X, complex<double>* Y, double delta, int isign)
{
        int fft_N = 2*N + 1;
        double* data = new double[fft_N];
        data[0] = 0.0;
        
        for (int i = 1; i < fft_N; i+=2)
        {  int m = (i-1)/2+1;
           data[i] = real(X[m-1]);
           data[i] *= real( exp( (double) isign * ii * pi * (double) m / (double) N) ); 
           data[i+1] = 0.0; //imag(X[m-1]);
        } 

        //--------------------------------------------//
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n = N << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
        //------------------------------//

        for(int i=1; i<fft_N; i+=2)
        { int n = (i-1)/2+1;
          Y[n-1] =  delta * complex<double>(data[i], data[i+1]);  
          Y[n-1] *= exp( (double) isign * ii * pi * (n-1.5) / (double) N );
        } 
        delete [] data;
}
#undef SWAP
//----------------------------------------------------------------.//

