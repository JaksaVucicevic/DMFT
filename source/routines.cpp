#include <cstdio>
#include "routines.h"
#include "GRID.h"
#include <vector>
#include <omp.h>
#include <iostream>
#include <cmath>
#include "nrutil.h"
// #include "mkl_types.h"
//#define MKL_Complex16 std::complex<double>
//#include "mkl.h"

#ifdef _OMP
#include <omp.h>
#endif

using namespace std;

double sign(double x)
{ //returns 1 if x positive, -1 if x negative, 0 if x=0
  if (x == 0.0) return 0;
  else return (x >= 0.0) ? +1.0 : -1.0;
}

int int_sign(double x)
{ if (x>=0.0) return 1; else return -1; };

double sqr(double x)
{ //returns square of x
  return x*x;
}

int pow(int base, int exp)
{
  int res = 1;
  for (int i=0; i<exp; i++)
    res *= base;
  return res;
}

complex<double> sqr(complex<double> x)
{
 return x*x;
}

/*double abs(double x)
{
  return (x>=0) ? x : -x;
}*/

/*double abs(complex<double> x)
{
  return  sqrt(   sqr( real(x) )
                + sqr( imag(x) ) );
}*/

void PrintFunc3D(const char* FileName, int N, complex<double>** Y, double* X)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i+=20)
    for (int j=0; j<N; j+=20)
      fprintf(f,"%.15le %.15le %.15le %.15le\n", X[i], X[j], real(Y[i][j]), imag(Y[i][j]));
  fclose(f);
}

void PrintFunc(const char* FileName, int N, int M, double** Y, double* X)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
  {
    fprintf(f, "%.15le", X[i]);
    for (int j=0; j<M; j++)
    {
       // loop through and store the numbers into the file
       fprintf(f, "%.15le", Y[i][j] );
    }
   fprintf(f, "\n");
  }
  fclose(f);
}

void PrintFunc(const char* FileName, int N, complex<double>* Y, double* X)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
    fprintf(f,"%.15le %.15le %.15le\n", X[i], real(Y[i]), imag(Y[i]));
  fclose(f);
}

void PrintFunc(const char* FileName, int N, complex<double>* Y)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
    fprintf(f,"%d %.15le %.15le\n", i, real(Y[i]), imag(Y[i]));
  fclose(f);
}

void PrintFunc(const char* FileName, int N, double* Y)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
    fprintf(f,"%d %.15le \n", i, Y[i]);
  fclose(f);
}


void PrintFunc(const char* FileName, int N, double* Y, double* X)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<N; i++)
    fprintf(f,"%.15le %.15le\n", X[i], Y[i]);
  fclose(f);
}

void PrintFunc(const char* FileName, std::vector< complex<double> > Y, std::vector<double> X)
{
  FILE *f;
  f = fopen(FileName, "w");
  for (int i=0; i<X.size(); i++)
    fprintf(f,"%.15le %.15le %.15le\n", X[i], real(Y[i]), imag(Y[i]));
  fclose(f);
}

void GetDimensions(const char* FileName, int &N, int &M)
{
 //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) perror ("Error opening file");

  //----count rows----//
  int i=0;
  char str[100000];
  while (!feof(f))
  {
    fgets ( str, 100000, f );
    i++;
  }
  N=i-1;

  //----count columns---//
  i=1;
  int j=0;
  while (str[i] != '\0')
  {  if ((str[i]!=' ')and(str[i+1]!='\0')and(str[i-1]==' ')) j++;
     i++;
  }
  M=j+1;

  //---close file-----//
  fclose(f);
}

void ReadFunc(const char* FileName, int &N, int &M, double** &X)
{
 GetDimensions(FileName, N, M);
 printf("N: %d M: %d\n", N, M);

 X = new double*[M];
 for (int i=0; i<M; i++)
   X[i] = new double[N];

 FILE *f;
 f = fopen(FileName, "r");

 for (int i=0; i<N; i++)
   for (int j=0; j<M; j++)
   { double Dummy;
     fscanf(f, "%le", &Dummy);
     X[j][i]=Dummy;
   }
 fclose(f);
 printf("uchito fajl");
}

void ReadFunc(const char* FileName, int &N, complex<double>* &Y, double* &X, bool PurelyReal)
{ // reads file formatted like 
  //  X  ReY(X) ImY(X)
  //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) perror ("Error opening file");
  
  //----count rows----//
  int i=0;
  char str[1000];  
  int prelines = 0;
  while (!feof(f))
  {  
    fgets ( str, 1000, f ); 
    if (str[0]=='\n') continue;
  
    string sline(str);   
    if ( ( sline.find("nan") != string::npos )
          or
         (str[0]=='#')
       )
       { prelines++; continue; };
    i++;
  }
  N=i-1;
  printf("N: %d, prelines: %d \n", N, prelines);
  fclose(f);
 
  X = new double[N];
  Y = new complex<double>[N];

  f = fopen(FileName, "r");

  for (int i=0; i<prelines; i++)
    fgets ( str, 1000, f ); 
   
  for (int i=0; i<N; i++)
  { double Dummy1,Dummy2,Dummy3;
    fscanf(f, "%le", &Dummy1);
    fscanf(f, "%le", &Dummy2);
    if (not PurelyReal) fscanf(f, "%le", &Dummy3);
    X[i]=Dummy1;
    Y[i]=complex<double>(Dummy2, (PurelyReal) ? 0.0 : Dummy3);
  }
  fclose(f);
}

void ReadFunc(const char* FileName, int &N, double* &Y, double* &X)
{ // reads file formatted like 
  //  X  ReY(X) ImY(X)
  //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) perror ("Error opening file");
  
  //----count rows----//
  int i=0;
  char str[1000];  
  int prelines = 0;
  while (!feof(f))
  {  
    fgets ( str, 1000, f ); 
    if (str[0]=='\n') continue;
  
    string sline(str);   
    if ( ( sline.find("nan") != string::npos )
          or
         (str[0]=='#')
       )
       { prelines++; continue; };
    i++;
  }
  N=i-1;
  printf("N: %d, prelines: %d \n", N, prelines);
  fclose(f);
 
  X = new double[N];
  Y = new double[N];

  f = fopen(FileName, "r");

  for (int i=0; i<prelines; i++)
    fgets ( str, 1000, f ); 
   
  for (int i=0; i<N; i++)
  { double Dummy1,Dummy2;
    fscanf(f, "%le", &Dummy1);
    fscanf(f, "%le", &Dummy2);

    X[i]=Dummy1;
    Y[i]=Dummy2;
  }
  fclose(f);
}

void ReadFunc(const char* FileName, int N, double* Y, double* X)
{ // reads file formatted like 
  //  X  ReY(X) ImY(X)
  //----open file---//
  FILE *f;
  f = fopen(FileName, "r");
  if (f==NULL) { printf("Error opening file"); return; }

  for (int i=0; i<N; i++)
  { double Dummy1,Dummy2;
    fscanf(f, "%le", &Dummy1);
    fscanf(f, "%le", &Dummy2);

    X[i]=Dummy1;
    Y[i]=Dummy2;
  }
  fclose(f);
}

bool FileExists(const char* FN){ FILE* f = fopen(FN,"r"); if (f==NULL) return false; else { fclose(f); return true; } };
//---------------- vectors and matrices--------------------//

void MultiplyByMatrix(int N, double* v, double** m)
{
  double* res = new double[N];
  for (int j=0; j<N; j++)
  { res[j] = 0;
    for (int i=0; i<N; i++)
      res[j] += m[i][j]*v[i];
  }
  for (int i=0; i<N; i++) v[i] = res[i];
  delete [] res;
}

void CreateRotationMatrix(double** m, int N, double angle, int* plane)
{
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      if (i==j) m[i][j] = 1.0;
      else m[i][j] = 0.0;
  m[plane[0]][plane[0]] = cos(angle);
  m[plane[1]][plane[1]] = cos(angle);
  m[plane[0]][plane[1]] = sin(angle);
  m[plane[1]][plane[0]] = -sin(angle);
}

void RotateVector(int N, double* v, double angle, int* plane)
{
  double** RotationMatrix = new double*[N];
  for (int i=0; i<N; i++)
    RotationMatrix[i] = new double[N];
  CreateRotationMatrix(RotationMatrix, N, angle, plane);
  MultiplyByMatrix(N, v, RotationMatrix);
  for (int i=0; i<N; i++) delete [] RotationMatrix[i];
  delete [] RotationMatrix;
}

void RotateVector2D(double* v, double angle)
{
  RotateVector(2, v, angle, (int []){0,1} );
}

#define NRANSI
#define TINY 1.0e-20

void lubksb(double **a, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}


void ludcmp(double **a, int n, int *indx, double *d)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv= new double[n+1];
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	delete [] vv;
}

void InvertMatrix(int N, double** A, double** invA, double &det)
{ 
  double** a = new double*[N+1];
  for(int i=0; i<N+1; i++) 
  { a[i] = new double[N+1];
    for(int j=0; j<N+1; j++) a[i][j] = ( (i>0) and (j>0) ) ? A[i-1][j-1] : 0.0; 
  }

  double *col = new double[N+1];
  int i,j,*indx;

  ludcmp(a,N,indx,&det); //Decompose the matrix just once.

  for(j=1;j<=N;j++) 
  { //Find inverse by columns.
    for(i=1;i<=N;i++) col[i]=0.0;
    col[j]=1.0;
    lubksb(a,N,indx,col);
    for(i=1;i<=N;i++) invA[i-1][j-1]=col[i];
  }
}

#undef TINY
#undef NRANSI










//------------------ integral routine ---------------------//

double TrapezIntegralMP(int N, double Y[], double X[])
{
#ifdef _OMP

  double sum = 0.0;
  double sum0 = Y[0]*(X[1]-X[0]) + Y[N-1]*(X[N-1]-X[N-2]);
  int Nt;
  double* psum = new double[8];
  #pragma omp parallel shared(psum)
  { Nt = omp_get_num_threads();
    int tid = omp_get_thread_num();
    psum[tid] = 0;
    for (int i=tid+1; i<N-1; i+=Nt)
      psum[tid]+=Y[i]*(X[i+1]-X[i-1]);
    //printf("inside proc %d, psum = %le\n",tid,psum[tid]);
  }
  for (int i = 0; i<Nt; i++)
  {  sum += psum[i];
     //printf("outside proc %d, psum = %le\n",i,psum[i]);
  }
  delete [] psum;
  return (sum+sum0)*0.5;

#else

  return TrapezIntegral(N,Y,X);

#endif
}

complex<double> TrapezIntegralMP(int N, complex<double> Y[], double X[])
{

#ifdef _OMP
  complex<double> sum0 = Y[0]*(X[1]-X[0]) + Y[N-1]*(X[N-1]-X[N-2]);

  complex<double> sum = 0.0;
  int Nt;
  complex<double>* psum = new complex<double>[8];
  #pragma omp parallel shared(psum)
  { Nt = omp_get_num_threads();
    int tid = omp_get_thread_num();
    psum[tid] = 0;
    for (int i=tid+1; i<N-1; i+=Nt)
      psum[tid]+=Y[i]*(X[i+1]-X[i-1]);
  }
  for (int i = 0; i<Nt; i++) sum += psum[i];
  delete [] psum;
  return (sum + sum0)*0.5;

#else

  return TrapezIntegral(N,Y,X);

#endif

}
//--------------------------- with sorted integrand -----------------------------------//
/*
double TrapezIntegral2D(int N, double** Y, double *X)
{
  for (int i=0; i<N; i++)
  for (int j=0; j<N; j++)
  { if (i==0) Y[i][j] *= 0.5*(X[1]-X[0]);
    if (i==N-1) Y[i][j] *= 0.5*(X[1]-X[0]);
    if ((i!=0)and(i!=N-1)) Y[i][j] *= 0.5 * (X[i+1]-X[i-1]);
    if (j==0) Y[i][j] *= 0.5*(X[1]-X[0]);
    if (j==N-1) Y[i][j] *= 0.5*(X[N-1]-X[N-2]);
    if ((j!=0)and(j!=N-1)) Y[i][j] *= 0.5 * (X[j+1]-X[j-1]);
  }
  int sqrN = N*N;
  double* integrand = new double[sqrN];
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      integrand[j*N+i] = Y[i][j];

  for (int i=0; i<sqrN-1; i++)
    for(int j=i+1; j<sqrN; j++)
      if (integrand[i]>integrand[j])
      {  double pom = integrand[j];
         integrand[j] = integrand[i];
         integrand[i] = pom;
      }
  double sum = 0.0;
  for (int i=0; i<sqrN; i++)
    sum+=integrand[i];
  return sum;
}

complex<double> TrapezIntegral2D(int N, complex<double>** Y, double *X)
{
  double** ReY = new double*[N];
  double** ImY = new double*[N];
  for (int i=0; i<N; i++)
  {
    ReY[i] = new double[N];
    ImY[i] = new double[N];
    for (int j=0; j<N; j++)
    {  ReY[i][j] = real(Y[i][j]);
       ImY[i][j] = imag(Y[i][j]);
    }
  }  
  double resum = TrapezIntegral2D(N, ReY, X);
  double imsum = TrapezIntegral2D(N, ImY, X);

  for (int i=0; i<N; i++)
  {  delete [] ReY[i];
     delete [] ImY[i];
  }
  delete [] ReY;
  delete [] ImY;  
  return complex<double>(resum,imsum);
}
*/
/*
double TrapezIntegral(int N, double Y[], double X[])
{

  for (int i=0; i<N; i++)
  { if (i==0) Y[i] *= 0.5*(X[1]-X[0]);
    if (i==N-1) Y[i] *= 0.5*(X[N-1]-X[N-2]);
    if ((i!=0)and(i!=N-1)) Y[i] *= 0.5 * (X[i+1]-X[i-1]);
  }
  for (int i=0; i<N-1; i++)
    for(int j=i+1; j<N; j++)
      if (abs(Y[i])>abs(Y[j]))
      {  double pom = Y[j];
         Y[j] = Y[i];
         Y[i] = pom;
      }
  double sum = 0.0;
  for (int i=0; i<N; i++)
    sum+=Y[i];

  return sum;
}


complex<double> TrapezIntegral(int N, complex<double> Y[], double X[])
{
  double* ReY = new double[N];
  double* ImY = new double[N];
  for (int i=0; i<N; i++)
  {
    ReY[i] = real(Y[i]);
    ImY[i] = imag(Y[i]);
  }  
  double resum = TrapezIntegral(N, ReY, X);
  double imsum = TrapezIntegral(N, ImY, X);
  delete [] ReY;
  delete [] ImY;  


  return complex<double>(resum,imsum);
}
*/
//--------------------------------------------------------------------------------------//


double TrapezIntegral(int N, double Y[], double X[])
{
  
  double sum = Y[0]*(X[1]-X[0]) + Y[N-1]*(X[N-1]-X[N-2]);
  for (int i=1; i<N-1; i++)
    sum+=Y[i]*(X[i+1]-X[i-1]);
  return sum*0.5;
}

complex<double> TrapezIntegral(int N, complex<double> Y[], double X[])
{
  complex<double> sum = Y[0]*complex<double>(X[1]-X[0]) +
                         Y[N-1]*complex<double>(X[N-1]-X[N-2]);
  for (int i=1; i<N-1; i++)
    sum+=Y[i]*complex<double>(X[i+1]-X[i-1]);

  return sum*0.5;
}


/*
double TrapezIntegral(std::vector<double> Y, std::vector<double> X)
{
  double sum = 0.0;
  int N = Y.size();
  sum += Y[0]*(X[1]-X[0]) + Y[N-1]*(X[N-1]-X[N-2]);
  for (int i=1; i<N-1; i++)
    sum+=Y[i]*(X[i+1]-X[i-1]);
  return sum*0.5;
}
*/
complex<double> TrapezIntegral(std::vector< complex<double> > Y, std::vector<double> X)
{
  complex<double> sum = 0.0;
  int N = Y.size();
  sum += Y[0]*complex<double>(X[1]-X[0]) +
         Y[N-1]*complex<double>(X[N-1]-X[N-2]);

  for (int i=1; i<N-1; i++)
    sum+=Y[i]*complex<double>(X[i+1]-X[i-1]);

  return sum*complex<double>(0.5);
}

double SmartIntegral(int N, double Y[], double X[], double x0, double Dx, double (*f)(double), int Nextra, double x_min, const char* FN)
{
  //                      +(xs-xc)Y[i1]                                   +(xc-xe)Y[i3]
  //    .    .    .         |#####|    .   .        .  |  .    .   .    |===| .    .    .
  //              i0        xc    xs  i1              x0          i2    xc xe i3   
  //                               <-------- Dx ------> <-------- Dx ------>
 
  double xs = x0-Dx;
  double xe = x0+Dx;

  if ( (xs > X[N-1]) or (xe<X[0]) ) return TrapezIntegral(N, Y, X);

  double sum1;
  if ( (xs>X[0])or(xe<X[N-1]) )
  { 
    int i1=N;
    int i2=N;
    int i0;
    int i3;

    if (xs<X[0]) i1=0;
    double* Ycopy = new double[N];
    for(int i=0; i<N; i++)
    { 
      if (i<N-1)
      { if ( (X[i]<=xs)and(X[i+1]>xs) )
        { i0=i;  i1=i+1; }
        if ( (X[i]<xe)and(X[i+1]>=xe) )
        { i2=i;  i3=i+1; }
      }

      if ((i>=i1)and(i<=i2)) 
      {  Ycopy[i]=0.0;
         //printf("i: %d, Ycopy=0\n",i);
      }
      else 
        Ycopy[i] = Y[i];
    }
    if (FN!=NULL)
    { char fullFN[300];
      sprintf(fullFN,"smartIntegral.regular.%s",FN);
      PrintFunc(fullFN,N,Ycopy,X);
    }
    
    double corr1, corr2;
    if(i1>0)
    {  double xc = 0.5*(X[i1]+X[i0]); 
       double w = xs - xc;
       corr1 = w * ((w<=0.0) ? Y[i0] : Y[i1]);
    }
    else
      corr1 = 0.0;

    if(i2<N)
    {  double xc = 0.5*(X[i2]+X[i3]); 
       double w = xc - xe;
       corr2 = w * ((w<=0.0) ? Y[i3] : Y[i2]);
    }
    else
      corr2 = 0.0;

    sum1 = TrapezIntegral(N, Ycopy, X) + corr1 + corr2;
    //printf("corr1: %.6f, corr2: %.6f, sum1: %.6f \n",corr1,corr2, sum1);
    delete [] Ycopy;

  }
  else sum1 = 0.0;


  double * Xextra = new double [Nextra];
  double * Yextra = new double [Nextra];
  
  GRID grid(Nextra, 0,  Dx+0.1, Dx, x_min);
  grid.assign_omega(Xextra);

  int i1=0;
  int i2=Nextra;
  int i0;
  int i3;

  for(int i=0; i<Nextra; i++)
  {  Xextra[i] += x0;
     if ( (Xextra[i]>=X[0]) and (Xextra[i]<=X[N-1]) )
       Yextra[i] = f(Xextra[i]); 
     else
       Yextra[i] = 0.0; 
  
     if (i<Nextra-1)
     { if ( (Xextra[i]<=X[0])and(Xextra[i+1]+x0>X[0]) )
       { i0=i;  i1=i+1; }
       if ( (Xextra[i]<X[N-1])and(Xextra[i+1]+x0>=X[N-1]) )
       { i2=i;  i3=i+1; 
         //printf("i2 = %d\n",i);
       }
     }
  }
  double corr1, corr2;
  if ((i1>0)and(Xextra[0]!=X[0]) )
  {  double xc = 0.5*(Xextra[i1]+Xextra[i0]); 
     double w = xc - X[0];
     corr1 = w * ((w<=0) ? Yextra[i1] : f(Xextra[i0]));
  }
  else
    corr1 = 0.0;

  if ( (i2<Nextra) and (Xextra[Nextra-1] != X[N-1]) )
  {  double xc = 0.5*(Xextra[i2]+Xextra[i3]); 
     double w = X[N-1] - xc;
     corr2 = w * ((w<=0) ? Yextra[i2] : f(Xextra[i3]));
  }
  else
    corr2 = 0.0;

  double sum2 = TrapezIntegral(Nextra, Yextra, Xextra) + corr1 + corr2;
  //printf("corr1: %.6f, corr2: %.6f, sum2: %.6f \n",corr1,corr2, sum2);
  if (FN!=NULL)
  { char fullFN[300];
    sprintf(fullFN,"smartIntegral.extra.%s",FN);
    PrintFunc(fullFN,Nextra,Yextra,Xextra);
  }
  
  delete [] Xextra;
  delete [] Yextra; 

  return sum1 + sum2; 
}





double EllipticIntegralFirstKind(double x)
{ //integration goes from 0 to 1
  //split integration interval
  int Nlog=1000;
  int Nlin=1000;
  int N = Nlog+Nlin;
  double* omega = new double[N];
  GRID grid(Nlog,Nlin,1.0,1.0e-1,1.0e-6);
  grid.assign_omega(omega);

  double* integrand = new double[N];
  for (int i=N/2; i<N; i++) integrand[i] = 0;
  for (int i=0; i<N/2; i++) integrand[i] = 1.0 / sqrt( (1-sqr(omega[i]+1.0))*(1-sqr((omega[i]+1.0)*x))  );

  double res = TrapezIntegral(N,integrand,omega);
  delete [] integrand;
  delete [] omega;
  return res;
}

complex<double> EllipticIntegralFirstKind(complex<double> x)
{ //integration goes from 0 to 1
  //split integration interval
  int Nlog=1000;
  int Nlin=1000;
  int N = Nlog+Nlin;
  double* omega = new double[N];
  GRID grid(Nlog,Nlin,1.0,1.0e-1,1.0e-14);
  grid.assign_omega(omega);

  complex<double>* integrand = new complex<double>[N];
  for (int i=N/2; i<N; i++) integrand[i] = 0;
  for (int i=0; i<N/2; i++) integrand[i] = 1.0 / (sqrt( complex<double> ( (1.0 - sqr(omega[i]+1.0))
                                                                              * ( 1.0 - sqr((omega[i]+1.0)*x) )
                                                                         )
                                                       )
                                                  );
//  /*char FN[50];
//  sprintf(FN,"elliptic.x%.3f",real(x));
//  PrintFunc(FN,N, integrand, omega);
  complex<double> res = TrapezIntegral(N,integrand,omega);
  delete [] integrand;
  delete [] omega;
  grid.~GRID();
  return res;
}

//------Sine Integral, from Num Recipes -----//

double SI(double x)
{
  //--------consts----------//
  double EPS = 3.0e-16;
  double EULER = 0.577215664901533;
  int MAXIT = 100;
  double PIBY2 = 1.570796326794897;
  double FPMIN = 1.0e-30;
  double TMIN = 2.0;
  complex<double> ONE =  complex<double>(1.0,0.0);
  //--------vars-------------//
  double si;
  int i,k,odd;
  double a,err,fact,sign,sum,sumc,sums,t,term;
  complex<double> h,b,c,d,del;
  t=fabs(x);
                                           // Special case.
  if (t == 0.0) {
      si=0.0;
//      *ci = -1.0/FPMIN;
//      return si;
  }
                                            //Evaluate continued fraction by modified
  if (t > TMIN) {
                                               // Lentz’s method (§5.2).
      b=complex<double>(1.0,t);
      c=complex<double>(1.0/FPMIN,0.0);
      d=h=ONE/b;
      for (i=2;i<=MAXIT;i++) {
          a = -(i-1)*(i-1);
          b+=complex<double>(2.0,0.0);
                                            //Denominators cannot be zero.
          d=ONE/( complex<double>(a,0)*d +b );
          c=b+complex<double>(a,0.0)/c;
          del=c*d;
          h=h*del;
          if (fabs(real(del)-1.0)+fabs(imag(del)) < EPS) break;
      }
      if (i > MAXIT) printf("cf failed in cisi\n");
      h=complex<double>(cos(t),-sin(t))*h;
//      *ci = -real(h);
      si=PIBY2+imag(h);
                                            //Evaluate both series simultaneously.
  } else {
                                            //Special case: avoid failure of convergence
      if (t < sqrt(FPMIN)) {
                                               // test because of underflow.
          sumc=0.0;
          sums=t;
      } else {
          sum=sums=sumc=0.0;
          sign=fact=1.0;
          odd=1;
          for (k=1;k<=MAXIT;k++) {
              fact *= t/k;
              term=fact/k;
              sum += sign*term;
              err=term/fabs(sum);
              if (odd) {
                  sign = -sign;
                  sums=sum;
                  sum=sumc;
              } else {
                  sumc=sum;
                  sum=sums;
              }
              if (err < EPS) break;
              odd=!odd;
          }
          if (k > MAXIT) printf("maxits exceeded in cisi\n");
      }
      si=sums;
//      *ci=sumc+log(t)+EULER;
  }
  if (x < 0.0) si = -(si);
  return si;
}


//-------------------- DOSes and Init Deltas -----------------//

double DOS(int DOStype, double t, double om, double U)
{
  switch (DOStype)
  {
    case DOStypes::SemiCircle:
            if ( sqr(2.0*t)<sqr(om) )
              return 0.0;
            else
              return sqrt( sqr(2.0*t)-sqr(om) )
                     / ( 2.0*pi*sqr(t) );
            break;
    case DOStypes::Gaussian:
            {
              double d = 1/(2*t*sqrt(pi))*exp( -sqr( om )/( 4*sqr(t) ) );
              return (d>0.001) ? d : 0;
            }
            break;
    case DOStypes::Insulator:
            {
              if ( sqr(2.0*t)<sqr(abs(om)-U/2.0) )
                return 0.0;
              else
                return ( ( sqr(2.0*t)>sqr(om-U/2.0) ) ? sqrt( sqr(2.0*t)-sqr(om-U/2.0) ) / ( 4.0*pi*sqr(t) ) : 0.0 )
                       +
                       ( ( sqr(2.0*t)>sqr(om+U/2.0) ) ? sqrt( sqr(2.0*t)-sqr(om+U/2.0) ) / ( 4.0*pi*sqr(t) ) : 0.0 );
            } break;
    case DOStypes::SquareLattice:
           { if (om==0.0) return 0.0;
             if (abs(om)>4.0*t) return 0.0;
             // !!!!!!!!!!!!!!! not good. use smart integral for Kelliptic
             if (abs(om)>=3.98*t) return (2.0/(4.0*t*sqr(pi)))*( 0.5*pi + pi*( 1.0-sqr(om/(4.0*t)))/8.0 );
             if (abs(om)<t/2.5) return log( 16.0/(sqr(om/(4.0*t)))) / (sqr(pi)*4.0*t);
             else return 2.0 / (sqr(pi)*4.0*t)
                         * EllipticIntegralFirstKind( sqrt( 1.0-sqr(om/(4.0*t))) );
           }  break;
    case DOStypes::CubicLattice:
           {  if (abs(om)>=6*t) return 0.0;
              int N = 6000;
              double* omega = new double[N];
              GRID grid(0, N, 3.141, 0.0, 0.0);
              grid.assign_omega(omega);
              double* d = new double[N];
              for (int i=0; i<N/2; i++) d[i]=0;
              #pragma omp parallel for
              for (int i=N/2; i<N; i++)
              {  complex<double> f = 4.0 * t / ( om - 2 * t * cos(omega[i]) + complex<double>(0.0, 0.001) );
                 //if (f>1.0)
                   d[i] = imag( f * EllipticIntegralFirstKind( f ) );
                 //else d[i]=0;
                 //if (om<1.0) printf("omega %.3f: f=%.3f\n",om,f);
              }
              /*if (om<1.0)
              {  char FN[50];
                 sprintf(FN,"integrand.w%.3f",om);
                 PrintFunc(FN,N,d,omega);
              }*/

              double Res = -1.0/pi * TrapezIntegral(N,d,omega);
              delete [] d;
              delete [] omega;
              grid.~GRID();
              return 1.0 / (2 * sqr(pi) * t)
                     * Res;
           }  break;
    case DOStypes::Uniform:
           {
             if ((om<-t)or(om>t)) return 0.0;
             else return 1.0/(2.0*t);
           } break;
  
    //Add here different DOS types!!!
    default: { printf("DOStype %d not implemented!", DOStype); exit(1); }
  }
}


void WriteCubicDosToFile()
{
  int N = 3000;
  double* omega = new double[N];
  double* dos = new double[N];
  GRID grid(0, N, 3.5, 0.0, 0.0);
  grid.assign_omega(omega);
  for(int i=N/2; i<N; i++)
  { dos[i] = DOS(DOStypes::CubicLattice, 0.5, omega[i]);
    dos[N-1-i]=dos[i];
    printf("Write2File: omega %.3f: %.3f\n", omega[i], dos[i]);
  }
  PrintFunc("CubicLatticeDOS",N,dos,omega);
  delete [] dos;
  delete [] omega;
}

void ReadDosFromFile(const char* FN, int N, double* omega, double* DOS)
{
  int n,m;
  double** dos;
  ReadFunc(FN, n, m, dos);
  for (int i = 0; i<N; i++)
    DOS[i] = interpl(n, dos[1], dos[0], omega[i]);
  for (int i=0; i<m; i++)
    delete [] dos[i];
  delete [] dos;
}

double interpl(int N, double* Y, double* X, double x)
{
  if ((x < X[0])or(x > X[N-1])) return 0;
  else
  { int i;
    for (i = 0; i < N; i++)
      if (X[i]>x) break;
    return Y[i-1] + (Y[i]-Y[i-1]) / (X[i]-X[i-1])
                                  * (x-X[i-1]);
  }
}

complex<double> interpl(int N, complex<double>* Y, double* X, double x)
{
  double* ReY = new double[N];
  double* ImY = new double[N];
  for(int i=0; i<N; i++)
  {
    ReY[i] = real(Y[i]);
    ImY[i] = imag(Y[i]);
  }

  double ReYx = interpl(N, ReY, X, x);
  double ImYx = interpl(N, ImY, X, x);

  delete [] ReY;
  delete [] ImY;

  return complex<double>(ReYx,ImYx);
}

//------------------------------------------------------------------------------//

void get_G_from_DOS(int DOStype, double t, int N, double* omega, complex<double>* G, double eta, double U)
{
  //valid only at half filling. dont use too small eta
  complex<double>* g = new complex<double>[N];

  for(int i=0; i<N; i++)
  {
    for(int j=0; j<N; j++)
      g[j] = DOS(DOStype, t, omega[j], U) / (omega[i]-omega[j]+ii*eta);
    G[i] = TrapezIntegral(N, g, omega);
  }
  
  delete [] g;
}

void get_Sigma_from_G(double t, int N, double* omega, complex<double>* G, complex<double>* Sigma)
{
  //valid only at half filling for Bethe lattice
  for(int i=0; i<N; i++)
    Sigma[i] = omega[i]-1.0/G[i]-sqr(t)*G[i];
  
}

complex<double> LS_get_G(int DOStype, double t, complex<double> com)
{
  complex<double> root = sqrt ( sqr(com) - 4.0 * sqr(t) );
  double sgn = sign(imag(root));
  return ( com - sgn*root) / ( 2.0 * sqr(t) ) ;
}

void InitG(int DOStype, double t, int N, double* omega, complex<double>* G)
{
  complex<double>* comega = new  complex<double>[N];
  for (int i=0; i<N; i++) comega[i] = omega[i];
  InitG(DOStype, t, N, comega, G);
  delete [] comega;
}

void InitG(int DOStype, double t, int N, complex<double>* omega, complex<double>* G)
{
  switch (DOStype)
  {
     case DOStypes::SemiCircle:
       for (int i = 0; i<N; i++)
         G[i] = LS_get_G(DOStype, t, omega[i]) ;
       break;
     default: printf("routines::InitG: NOT IMPLEMENTED!!!");
  }
}



void InitDOS(int DOStype, double t, int N, double* omega, double* dos, double U)
{
  for (int i=0; i<N; i++) dos[i] = DOS(DOStype, t, omega[i], U);
}

void InitDelta(int DOStype,
               int N,
               double V,
               double mu,
               double eta,
               double t,
               double* omega,
               complex<double>* Delta,
               const char* FNDos)
{
  //printf("-- INFO -- routines: Making Delta: t=%.3f\n",t);
  // prepare dos

//------------------Get DOS----------------------//
  double* dos = new double[N];


// TODO : THIS NEEDS TO BE FIXED

//  if (FNDos=="")
// // // // // // // // // // // // // // // // // // // // // // // // // // // // //     #pragma omp parallel for
    for(int i=0; i<N; i++)
      dos[i] = DOS(DOStype, t, omega[i]);
//  else
//    ReadDosFromFile(FNDos, N, omega, dos);

  PrintFunc("InitDelta.DOS",N,dos,omega);

//-----------------------------------------------//
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // //   #pragma omp parallel for
  for (int i=1; i<N-1; i++)
  {
    //treat integrand carefully
    double D = 0.0;
    complex<double> LogTerm = 0.0;
    D = dos[i];

    LogTerm =  ((i==0)or(i==N-1)) ? 0.0
                                  : complex<double>(D, 0.0) * log( complex<double> ( (mu + omega[i] + omega[N-1])
                                                                                    /(mu + omega[i] - omega[N-1]) ) );
    //printf("LogTerm(%.3f) = %.3f, %.3f\n",omega[i], real(LogTerm),imag(LogTerm));
    //create integrand array
    complex<double>* d = new complex<double>[N];
    for (int j=0; j<N; j++)
    {  d[j] = (i==j) ? 0.0 : complex<double>(dos[j] - D, 0.0)
                            / complex<double>( mu + omega[i] - omega[j] );
       //if (i==j) printf("d(i=j) = %.3f, %.3f\n", real(d[j]), imag(d[j]) );
    }

    //integrate to get G
    Delta[i] = conj( sqr(V)*(TrapezIntegral(N, d,omega) + LogTerm) ) ;
    delete [] d;
  }
  Delta[0] = Delta[1];
  Delta[N-1] = Delta[N-2];
  PrintFunc("DeltaMade",N,Delta,omega);

//---------------------------------------------------------------------------//

  delete [] dos;
  //printf("-- INFO -- routines: Intgral Delta: %.6f\n", imag(TrapezIntegral(N,Delta,omega)));
}

void InitDeltaFromSIAMOutputFile(const char* FN, int N, double* omega, complex<double>* Delta)
{
  int n,m;
  double** output;
  ReadFunc(FN, n, m, output);
  for (int i = 0; i<N; i++)
    Delta[i] = complex<double>( interpl(n, output[2], output[0], omega[i]),
                              interpl(n, output[3], output[0], omega[i]) );
  for (int i=0; i<m; i++)
    delete [] output[i];
  delete [] output;
}

void InitInsulatorDelta(double U,
                        int N,
                        double t,
                        double* omega,
                        complex<double>* Delta)
{ printf("Making insulator Delta from Hubbar I approx: U=%.3f, t=%.3f\n",U,t);
  for (int i=0; i<N; i++)
  {
    double Gat = 0.5 * ( 1.0/(omega[i] + U/2.0) + 1.0/(omega[i] - U/2.0) );
    complex<double> Sqr = sqrt( complex<double>(1.0)/sqr(Gat)  - 4.0*sqr(t) );
    complex<double> G = (1.0/Gat + real(Sqr) < 0) ? (1.0/Gat + conj(Sqr) ) / (2*sqr(t))
                                                  : (1.0/Gat - Sqr ) / (2*sqr(t));
    Delta[i] = sqr(t)*G;
  }
  PrintFunc("initdelta",N,Delta,omega);
   printf("Intgral Delta: %.6f", imag(TrapezIntegral(N,Delta,omega)));
}

//-----------------------------------------------------------------------------//
/*
void InitArrays(int N, double*** a, const int* Ns)
{
  for (int i=0; i<N; i++)
    a[i][0] = new double[Ns[i]];
}

void InitArrays(int N, double*** a, const int* Ns, double (**func)(int))
{
  for (int i=0; i<N; i++)
  {   a[i][0] = new double[Ns[i]];
      if (func[i]!=NULL)
        for (int j=0; j<Ns[i]; j++)
        {  a[i][0][j]=func[i](j);
           printf("j: %d a: %f\n",j,a[i][0][j]);
        }
  }
}

void ReleaseMemory(int N, double*** a)
{
  for (int i=0; i<N; i++)
    delete [] a[i][0];

}


void Transpose(complex<double> *Transposed, complex<double> *A ,int n)
{
  int i,j;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
    {
	real(Transposed[i*n+j]) = real(A[i+j*n]);
	imag(Transposed[i*n+j]) = imag(A[i+j*n]);
    }
}

void TransposeS(const char UL, MKL_Complex16 *Transposed, complex<double> *A ,int n)
// Transposing regular Hermitian matrix to upper or lower triangular
{
  int i,j,ul;
  if ((UL == 'U') or (UL == 'u'))
    ul = 1;
  if ((UL == 'L') or (UL == 'l'))
    ul = -1;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if (i*ul >= j*ul )
      {
	Transposed[i*n+j].real = real(A[i+j*n]);
	Transposed[i*n+j].imag = imag(A[i+j*n]);
      }
      else
      {
	Transposed[i*n+j].real = 0.0;
	Transposed[i*n+j].imag = 0.0;
      }
}
void TransposeSB(const char UL, complex<double> *Transposed, MKL_Complex16 *A ,int n)
// Transposing Hermitian upper or lower triangular matrix to regular
{
  int i,j,ul;
  if ((UL == 'U') or (UL == 'u'))
    ul = 1;
  if ((UL == 'L') or (UL == 'l'))
    ul = -1;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      if (i*ul <= j*ul )
      {
	real(Transposed[i*n+j]) = A[i+j*n].real;
	imag(Transposed[i*n+j]) = A[i+j*n].imag;
	real(Transposed[i+n*j]) = A[i+j*n].real;
	imag(Transposed[i+n*j]) = A[i+j*n].imag;
      }
}

void Invert_Matrix(complex<double>* invA, complex<double>* A, int N)
{
  int LWORK=4*N*N;
  int *permutations;
  MKL_Complex16 *WORK, *tempA;
  tempA = new MKL_Complex16 [N*N];
  permutations = new int[2*N];
  WORK = new MKL_Complex16 [4*N*N];

  int INFO=0;
  char UL = 'L';
  TransposeS(UL,tempA,A,N);
  //cout << tempA[0].real << "  " << tempA[0].imag << endl;
  int *IPIV;
  IPIV = new int[N];
  //zhetrf(&UL, &N, tempA, &N, IPIV, WORK, &LWORK, &INFO);
  zsytrf_(&UL, &N, tempA, &N, IPIV, WORK, &LWORK, &INFO);
  //zgetrf_( &N, &N, tempA , &N, IPIV, &INFO );
  if (INFO != 0)
  {
    cout << "ComplexMatrixInverse: Error at zhetrf INFO = " << INFO; exit(0);
  }
  //zhetri(&UL, &N, tempA, &N, IPIV, WORK, &INFO);
  //zsytri2_(&UL, &N, tempA, &N, IPIV, WORK, &LWORK, &INFO);
  zsytri_(&UL, &N, tempA, &N, IPIV, WORK, &INFO);
  //zgetri_( &N, tempA , &N, IPIV, WORK, &LWORK, &INFO );
  if (INFO != 0)
  {
    cout << "ComplexMatrixInverse: Error at zhetri  \n"; exit(0);
  }
  TransposeSB(UL,invA,tempA,N);

  delete [] WORK;
  delete [] tempA;
  delete [] permutations;
  delete [] IPIV;

}
*/
