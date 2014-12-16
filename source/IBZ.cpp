#include "IBZ.h"
#include "../source/routines.h"
#include "../source/arrayInitializers.h"
#include <omp.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

IBZ::IBZ(int IBZtype, int Nx, int Ny, double a)
{
  t = 0.25;
  this->a = a;
  this->IBZtype = IBZtype;
  this->Nx = Nx;
  this->Ny = Ny;
  Ntotal = Nx*Ny;
  Initialize();
}
IBZ::~IBZ()
{
  ReleaseMemory();
}

void IBZ::Initialize()
{
  kx=new double[Nx];
  ky=new double[Ny];
  for(int i=0; i<Nx; i++)
    kx[i] = get_kx(i);
  for(int j=0; j<Ny; j++)
    ky[j] = get_ky(j);

  summand = Array2D< complex<double> >(Nx, Ny);
  epsilon = Array2D<double>(Nx, Ny);
  velocity = Array2D <double>(Nx, Ny);
  for(int i=0; i<Nx; i++)
  for(int j=0; j<Ny; j++)
  {
    velocity[i][j] = get_velocity(kx[i],ky[j]);
    epsilon[i][j] = get_epsilon(kx[i],ky[j]);
  }
  
}

void IBZ::ReleaseMemory()
{
  FreeArray2D< complex<double> >(summand, Nx);
  FreeArray2D<double>(epsilon, Nx);
  FreeArray2D<double>(velocity, Nx);
}

void IBZ::SetLatticeSpacing(double a)
{
  ReleaseMemory();
  this->a = a;
  Initialize();
}  

complex<double> IBZ::sum()
{
  switch (IBZtype) 
  {
    case IBZtypes::SquareLattice: return 0;
    case IBZtypes::TriangularLattice: 
      #ifdef _OMP
       int Nt;      //printf("Nt: %d\n",Nt);
       complex<double> res = 0; 
       complex<double>* pres = new complex<double>[16];
       #pragma omp parallel shared(pres)
       { Nt = omp_get_num_threads(); 
         //printf("Nt: %d\n",Nt);
         int tid = omp_get_thread_num();
         //printf("tid: %d\n",tid);
         pres[tid] = 0;
         //for (int i=tid; i<Nx; i+=Nt)
         
         for (int j=tid; j<Ny; j+=Nt)
         for(int i=0; i<Nx; i++)
           pres[tid]+=summand[i][j];
       }
       for(int l=0; l<Nt; l++) res+=pres[l];
       delete [] pres;
      #else
       complex<double> res = 0; 
       for(int i=0; i<Nx; i++)
       for(int j=0; j<Ny; j++)
         res+=summand[i][j];
      #endif
      return res/(double) Ntotal; 
    default: return 0;
  }
}
void IBZ::PrintToFile(const char* FN)
{
  FILE* f = fopen(FN,"w");
  for(int i=0; i<Nx; i++)
  { for(int j=0; j<Ny; j++)
      fprintf(f,"%.15le %.15le %.15le %.15le %.15le\n", kx[i], ky[j], epsilon[i][j], real(summand[i][j]), imag(summand[i][j]));
    fprintf(f,"\n");
  } 
  fclose(f);
}

double IBZ::get_kx(int i)
{
  switch (IBZtype) 
  {
    case IBZtypes::SquareLattice: return 0;
    case IBZtypes::TriangularLattice: 
      return i * 4.0 * pi / ( Nx * a ); 
    default: return 0;
  }
}
double IBZ::get_ky(int j)
{
  switch (IBZtype) 
  {
    case IBZtypes::SquareLattice: return 0;
    case IBZtypes::TriangularLattice: 
      return    j * 2.0 * pi
             / ( sqrt(3.0) * Nx * a ) ; 
    default: return 0;
  }
}

double IBZ::get_epsilon(double kx, double ky)
{
  switch (IBZtype) 
  {
    case IBZtypes::SquareLattice: return 0;
    case IBZtypes::TriangularLattice: 
      return -2.0 * t * ( cos(kx*a) + 2.0 * cos(kx*a/2.0) * cos(sqrt(3.0)*ky*a/2.0) ); 
    default: return 0;
  }
}

double IBZ::get_velocity(double kx, double ky)
{
  switch (IBZtype) 
  {
    case IBZtypes::SquareLattice: return 0;
    case IBZtypes::TriangularLattice: 
      return 2.0 * t * a * ( sin(kx*a) + sin(kx*a/2.0) * cos(sqrt(3.0)*ky*a/2.0) ); 
    default: return 0;
  }
}


