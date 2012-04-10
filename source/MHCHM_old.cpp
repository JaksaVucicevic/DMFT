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
 
  rhot = new complex<double>[N];
  Sigmat = new complex<double>[N]; 

  get_G0();
  r->PrintResult("after_G0");
  printf("got G0...\n");
  
  get_rhot();
  printf("got rhot...\n");

  get_Sigmat();
  printf("got Sigmat...\n");

  get_Sigma();
  printf("got Sigma...\n");

  get_G();
  printf("got G...\n");

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
  int fft_N = 2*N + 1;
  double* fft_rho = new double[fft_N];
  double domega = grid->get_domega();
  fft_rho[0]=0.0;
  for(int i=1; i<fft_N; i+=2)
  {  
    int n = (i-1)/2+1;
    fft_rho[i] = - (1.0 / pi) * imag(r->G0[n-1]); 
    fft_rho[i] *=  */ real( exp(-ii * pi * (double) n / (double) N) );
    fft_rho[i+1] = 0.0; 
  }  

  PrintFunc("fft_rho_before", fft_N, fft_rho);

  FastFourier(fft_rho, N, -1);

  for(int i=1; i<fft_N; i+=2)
  { int m = (i-1)/2+1;
    rhot[m-1] = domega * complex<double>(fft_rho[i], fft_rho[i+1]);
    rhot[m-1] *= exp( - ii * pi * (m-1.5) / (double) N );
  }

  PrintFunc("fft_rho_after", fft_N, fft_rho);
  PrintFunc("rhot", N, rhot);

  delete [] fft_rho;

}

void MHCHM::get_Sigmat()
{
  printf("U: %f", U);
  for (int j = 0; j<N/2; j++)
    Sigmat[j] = pow(U, 2.0) * pow(rhot[j], 3.0);
  for (int j = N/2; j<N; j++)
    Sigmat[j] = 0.0;
  PrintFunc("Sigmat", N, Sigmat);
}

void MHCHM::get_Sigma()
{
  int fft_N = 2*N + 1 ;
  double* fft_sigma = new double[fft_N];
  double dt = 2.0 * pi / (N*grid->get_domega());

  fft_sigma[0] = 0.0;
  for (int i = 1; i < fft_N; i+=2)
  {  int m = (i-1)/2+1;
     fft_sigma[i] = real(Sigmat[m-1]);
     fft_sigma[i] *= real( exp(ii * pi * (double) m / (double) N) ); 
     fft_sigma[i+1] = 0.0;
  }

  PrintFunc("fft_sigma_before", fft_N, fft_sigma);
  FastFourier(fft_sigma, N, 1);
  
  for(int i=1; i<fft_N; i+=2)
  { int n = (i-1)/2+1;
    r->Sigma[n-1] = dt * complex<double>(fft_sigma[i], fft_sigma[i+1]);  
    r->Sigma[n-1] *= exp( ii * pi * (n-1.5) / (double) N );
    r->Sigma[n-1] *= - 2.0 * ii /** (double) N*/ ;
    if (imag(r->Sigma[n-1]) > 0.0) 
      r->Sigma[n-1] = complex<double>(real(r->Sigma[n-1]),0.0);
    if (n-1 >= N/2 ) r->Sigma[n-1] = 0.0;
  }

  PrintFunc("fft_sigma_after", fft_N, fft_sigma);
  //grid->KramarsKronig(r->Sigma);
  r->PrintResult("after_Sigma");
  
  delete [] fft_sigma;
}

void MHCHM::get_G()
{
  /*for (int i = 0; i<N; i++)
    r->G[i] = 1.0 / (r->omega[i] - r->Sigma[i] - r->Delta[i]);*/

  for (int i = 0; i<N; i++)
  { complex<double> zom = r->omega[i] - r->Sigma[i] + ii*SIAMeta;
    complex<double> root = sqrt ( sqr(zom) - 4.0 * sqr(t) );
    double sgn = sign(imag(root));
    r->G[i] = (zom - sgn*root) / ( 2.0 * sqr(t) ) ;
    if (i>N/4) r->G[i] = 0.0;
  }
}

void MHCHM::assign_NIG(Result* r)
{
  int N = r->grid->get_N();
  for (int i = 0; i<N; i++)
  { complex<double> zom = r->omega[i] + ii*SIAMeta;
    complex<double> root = sqrt ( sqr(zom) - 4.0 * sqr(t) );
    double sgn = sign(imag(root));
    r->Delta[i] = sqr(t) * (zom - sgn*root) / ( 2.0 * sqr(t) ) ;
  }
}

//----------Fast Fourier-----------------------------//
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void MHCHM::FastFourier(double data[], int N, int isign)
{
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
}
#undef SWAP
//----------------------------------------------------------------.//

