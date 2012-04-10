#include <cmath>
#include <complex>
#include <iostream>
#include "routines.h"
#include "Result.h"
#include "GRID.h"
#include "SIAM.h"
#include "TMT.h"
#include "Input.h"

#ifdef _MPI
#include "mpi.h"
#endif

#ifdef _OMP
#include <omp.h>
#endif

void TMT::Defaults()
{
  W = 0.1;
  Distribution = Distributions::Uniform;
  Nimp = 10;
  CHM::UseBethe = true;

  AverageNt = 8;
  SiamNt = 1;
  KramarsKronigNt = 8;

  ExitSignal = -1000.0;
}

TMT::TMT() : CHM()
{
  Defaults();
  mu0grid = new double[Nimp];
}

TMT::TMT(const char* ParamsFN) : CHM(ParamsFN)
{
  Defaults();

  this->ParamsFN.assign(ParamsFN);
  cout << "-- INFO -- TMT: Params File Name set to:" << this->ParamsFN << endl;

  Input input("params");
  
  input.ReadParam(W,"TMT::W");
  input.ReadParam(Distribution,"TMT::Distribution");
  input.ReadParam(Nimp,"TMT::Nimp");
  input.ReadParam(AverageNt,"TMT::AverageNt");
  input.ReadParam(SiamNt,"TMT::SiamNt");
  input.ReadParam(KramarsKronigNt,"TMT::KramarsKronigNt");
  
  mu0grid = new double[Nimp];
  for (int i = 0; i<Nimp; i++) mu0grid[i] = 0;
  //----------------------------------------------------------------------------------------------------//
  // NEVER FORGET TO SET VALUES TO A NEWLY INITIALIZED DYNAMIC ARRAY. BY DEFAULT THESE VALUES ARE NOT 0 //
  //----------------------------------------------------------------------------------------------------//
}

TMT::~TMT() 
{
  siam->~SIAM();
  delete [] mu0grid;
}

void TMT::SetWDN(double W, int Distribution, int Nimp)
{
  this->W = W;
  this->Distribution = Distribution;
  this->Nimp = Nimp;
  delete [] mu0grid;
  mu0grid = new double[Nimp];
  for (int i = 0; i<Nimp; i++) mu0grid[i] = 0;
}

void TMT::SetUseBethe(bool UseBethe)
{
  printf("-- INFO -- TMT: Must use Bethe lattice in TMT\n");

}

void TMT::SendExitSignal()
{
#ifdef _MPI
  int Nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &Nproc);
  for (int p = 1; p<Nproc; p++)
     MPI_Send(&ExitSignal, 1, MPI_DOUBLE, p, 99, MPI_COMM_WORLD);
#endif
}

double TMT::P(double epsilon)
{
  switch (Distribution)
  {
    case Distributions::Uniform : return 1.0/W; 
    case Distributions::Gaussian :
    { return 0; //TODO 
    } 
    default : return 0;
  }
}

void TMT::MakeEgrid()
{
  Egrid = new double[Nimp];
  for (int i=0; i<Nimp; i++)
    Egrid[i] = (W!=0.0) ? - W/2.0 + i * W / ( Nimp - 1.0 ) : 0;
}

void TMT::Avarage(Result** R)
{
#ifdef _OMP
  omp_set_num_threads(AverageNt);
#endif

  #pragma omp parallel shared(R)
  { 
    #pragma omp parallel for 
    for (int i=0; i<N; i++)
    { 
      double dosi = 0;
      for (int j=0; j<Nimp; j++)
        dosi += (1.0/Nimp) * log(R[j]->DOS[i]);        
      r->DOS[i] = exp(dosi);
      r->G[i] = complex<double>(0.0, -pi * r->DOS[i]); 

      double dosmedi = 0; 
      for (int j=0; j<Nimp; j++)
        dosmedi += R[j]->DOS[i];        
      r->DOSmed[i] = dosmedi/Nimp;    
    }
  } //end of parallel

#ifdef _OMP
  omp_set_num_threads(KramarsKronigNt);
#endif

  r->grid->KramarsKronig(r->G);
}

void TMT::PrepareResult(Result* R, double mu, double mu0, double* ReDelta, double* ImDelta)
{
  grid->assign_omega(R->omega);
      
  R->mu = mu;
  R->mu0 = mu0;
    
  for (int j=0; j<N; j++)
    R->Delta[j] = complex<double>(ReDelta[j],ImDelta[j]);
}

bool TMT::DoSIAM(Result* R, double epsilon)
{
  //SIAM siam(ParamsFN.c_str());
  //SIAM siam;
  siam->SetUTepsilon(U, T, epsilon);
  siam->SetIsBethe(UseBethe);
  siam->SetBroadening(SIAMeta); 
  //siam.SetBroydenParameters(50, 1e-12);

  return siam->Run(R);
}


void TMT::Slave(int myrank)
{
#ifdef _MPI  
  MPI_Status status;
  while (true)
  { 
    double mu; 
    int Nepsilons;
    MPI_Recv(&mu, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
    if (mu == -1000.0) break;
    MPI_Recv(&Nepsilons, 1, MPI_INT, 0, 99, MPI_COMM_WORLD, &status);

    printf("================== PROC %d =================\n",myrank);

    double* epsilons = new double[Nepsilons];
    double* mu0s = new double[Nepsilons];
    
    MPI_Recv(epsilons, Nepsilons, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
    MPI_Recv(mu0s, Nepsilons, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);

    double* ReDelta = new double[N];
    double* ImDelta = new double[N];

    MPI_Recv(ReDelta, N, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
    MPI_Recv(ImDelta, N, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
 
    MPI_Recv(&U, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
    MPI_Recv(&T, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
    
    Result** R = new Result*[Nepsilons];  

    bool Error = false;
    for (int i=0; i<Nepsilons; i++)
    {  
       R[i] = new Result(grid);
       PrepareResult(R[i], mu, mu0s[i], ReDelta, ImDelta);

       omp_set_num_threads(SiamNt);
       printf("PROC %d ::: ",myrank);
       //char FN[50];
       //sprintf(FN,"recieved.eps%.3f",epsilons[i]);
       //R[i]->PrintResult(FN);
       if (!Error)
         Error = DoSIAM(R[i], epsilons[i]);
    }

    for (int i=0; i<Nepsilons; i++)
    {  if (Error)
       {  double ErrorSignal = -1000.0;
          MPI_Send(&ErrorSignal, 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
       }
       else
         MPI_Send(&(R[i]->mu0), 1, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
       MPI_Send(R[i]->DOS, N, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
       R[i]->~Result();
    }

    delete [] R;
    delete [] epsilons;
    delete [] mu0s;
    delete [] ReDelta;
    delete [] ImDelta;

  }
#endif
}
//========================= MPI ===============================//
bool TMT::SolveSIAM()
{
/*  #ifdef _MPI
  printf("---- PARALLEL MPI\n");
  #endif
  
  #ifdef _OMP
  printf("---- PARALLEL OMP\n");
  #endif
*/

#ifdef _MPI 

  bool Error = false;

  MPI_Status status;
  
  MakeEgrid();

  //r->PrintResult("initial");

  double* ReDelta = new double[N];
  double* ImDelta = new double[N]; 
  for (int i=0; i<N; i++)
  { ReDelta[i] = real(r->Delta[i]);
    ImDelta[i] = imag(r->Delta[i]);
  } 

  int Nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &Nproc);

  printf("==== MPI ==== NUMBER OF PROCESSES: %d",Nproc); 
  
  int Nsp = Nimp % Nproc;
  int Nbare = Nimp / Nproc;  
  int* Nepsilons = new int[Nproc];
  double** epsilons = new double*[Nproc];
  double** mu0s = new double*[Nproc];
  
  for (int p=0; p<Nproc; p++)
  { 
    Nepsilons[p] = Nbare;
    if (p<=Nsp-1) Nepsilons[p]++;
    epsilons[p] = new double[Nepsilons[p]];
    mu0s[p] = new double[Nepsilons[p]];
    for (int i=0; i<Nepsilons[p]; i++)
    {  epsilons[p][i] = Egrid[p+i*Nproc];
       mu0s[p][i] = mu0grid[p+i*Nproc];
    }

    if (p!=0)
    {  MPI_Send(&(r->mu), 1, MPI_DOUBLE, p, 99, MPI_COMM_WORLD); 
       MPI_Send(&(Nepsilons[p]), 1, MPI_INT, p, 99, MPI_COMM_WORLD); 
       MPI_Send(epsilons[p], Nepsilons[p], MPI_DOUBLE, p, 99, MPI_COMM_WORLD);
       MPI_Send(mu0s[p], Nepsilons[p], MPI_DOUBLE, p, 99, MPI_COMM_WORLD);
       MPI_Send(ReDelta, N, MPI_DOUBLE, p, 99, MPI_COMM_WORLD);  
       MPI_Send(ImDelta, N, MPI_DOUBLE, p, 99, MPI_COMM_WORLD);  
       MPI_Send(&U, 1, MPI_DOUBLE, p, 99, MPI_COMM_WORLD);  
       MPI_Send(&T, 1, MPI_DOUBLE, p, 99, MPI_COMM_WORLD);  
    }

  }

  //----------- DO CALC and collect data ------------//
  Result** R = new Result*[Nimp];
  printf("============== === === === ==== MASTER: COLLECTING DATA...\n");
  for(int p=0; p<Nproc; p++)
  {  //printf("---------- Proc %d\n",p);
     for(int i=0; i<Nepsilons[p]; i++)
     { 
       int imp = p+i*Nproc;
       //printf("-------------- Imp %d\n",imp); 
       R[imp] = new Result(*r);
       if (p==0)
       { //PrepareResult(R[imp], r->mu, mu0s[p][i], ReDelta, ImDelta); 
         R[imp]->mu0 = mu0s[p][i];         
         omp_set_num_threads(SiamNt);
         printf("PROC 0 ::: ");
         if (!Error)
           Error = DoSIAM(R[imp], epsilons[p][i]);
         mu0grid[imp] = R[imp]->mu0;
       } 
       else
       {
         MPI_Recv(&(mu0grid[imp]), 1, MPI_DOUBLE, p, 99, MPI_COMM_WORLD, &status);
         if (!Error) Error = (mu0grid[imp] == -1000.0);

         MPI_Recv(R[imp]->DOS, N, MPI_DOUBLE, p, 99, MPI_COMM_WORLD, &status);
       }
       
       /*char IFN[50]; 
       sprintf(IFN,"collected.%d",imp);
       R[imp]->PrintResult(IFN);*/
     }
  }
  printf("DONE!!!");

  //---- release memory----//
  for (int p=0; p<Nproc; p++)
  { delete [] epsilons[p];
    delete [] mu0s[p];
  }
  delete [] epsilons;
  delete [] mu0s;
  delete [] ReDelta;
  delete [] ImDelta;
  delete [] Nepsilons;

  //------- average results---------//
  Avarage(R);
  r->PrintResult("Averaged");


  //----- release memory-------//
  for (int i=0; i<Nimp; i++)
   R[i]->~Result();
  
  delete [] R;
  delete [] Egrid;

  return Error;

#else

  bool Error = false;

  MakeEgrid();

  Result** R = new Result*[Nimp];
 
  for (int i=0; i<Nimp; i++)
  {
     R[i] = new Result(*r);
     R[i]->mu0 = mu0grid[i];

     Error  = DoSIAM(R[i], Egrid[i]);
        
     mu0grid[i] = R[i]->mu0;  
  }

  Avarage(R);

  // Release Memory
  for (int i=0; i<Nimp; i++)
   R[i]->~Result(); 
  delete [] R;

  delete [] Egrid;

  return Error;

#endif
}



