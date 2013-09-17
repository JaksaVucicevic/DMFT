//*************************************************//
//           Modified Broyden Solver for           //
//             systems of N equations              //
//                                  by JV 2010     //
//*************************************************//

#include <cstdio>
#include <complex>
#include <cstdlib>

using namespace std;

/*namespace BroydenGlobal
{
  bool PrintOutIterations = true;
};*/

class Broyden
{
  private:
    bool Initialized;   

    //--user set parameters--//
    int N;			//dimesion of the problem (determines the size of V, Vold, F, Fold, DV, DF and U).
    int MAX_ITS;		//maximum number of iterations (determines the size od c, U and Beta). Iterations start from 1.
    double alpha;		//mixing parameter
    double omega0;
    double Accr;

    //--LastReset--//
    int LastReset;
    //--storage arrays--//
    complex<double>* c;		//MAX_ITS x 1
    complex<double>** U;		//MAX_ITS x N
    complex<double>** Beta;	//MAX_ITS x MAX_ITS
    complex<double>** A;		//MAX_ITS x MAX_ITS
    complex<double>* V;		// N x 1
    complex<double>* Vold;	// N x 1
    complex<double>* F;		// N x 1
    complex<double>* Fold;	// N x 1
    complex<double>** DV;	// MAX_ITS x N
    complex<double>** DF;	// MAX_ITS x N

    void add_Ds(int it);
    void add_As(int it);
    bool get_Betas(int it);
    void get_Us(int it);
    void get_cs(int it);
    complex<double> CorrTerm(int it, int i);

    //--Algebra--//
    int KDelta(int i, int j); 							//returns Kronecker delta (i,j)
    complex<double> Multiply(complex<double> A[], complex<double> B[], int n); 	//returns Dot product of 2 complex (0..n-1) vectors
    bool ludcmp(complex<double> **a, int n, int *indx);				//LU decomposes **a
    void lubksb(complex<double> **a, int n, int *indx, complex<double> b[]);	//backtrace method
    bool InverseMatrix(complex<double> **a, complex<double> **y, int n);	//returns inverse of a (1..n x 1..n) in y  

    void PrepareArrays();
    void ReleaseMemory();
  public:
    //--user interface--//
    Broyden();
    void SetParameters(int N, int MAX_ITS, double alpha, double omega0, double Accr);
    void TurnOn(int it);
    void TurnOff();
    void Reset(int it);
    int CalculateNew(complex<double> Vnew[], int it);
    double CurrentDiff;
};


/***********************************************************************
  Use this function to solve a system of N equations formulated like:
                       \vec{F}(\vec{V})=\vec{V}

  It should be implemented as a public member Function of a class, e.g:

  class T { public: void F(complex<double>* V); };

  N - number of equations (length of V)
  MAX_ITS - maximum number of broyden iterations
  Accr - requiered accuracy
  func - F
  obj - the instantianted object of which member function will be called
  V - initial guess on input, result on output 
************************************************************************/
template <class T> 
bool UseBroyden(int N, int MAX_ITS, double Accr, 
                void (T::*func)(complex<double>*), T* obj, complex<double>* V)
{
  bool b = false;

  //init broyden
  Broyden B;
  B.SetParameters(N, MAX_ITS, 0.99, 0.01, Accr);
  B.TurnOn(0);

  //save initial guess
  complex<double>* Vinit = new complex<double>[N];
  for (int i=0; i<N; i++) Vinit[i] = V[i];

  //------------- iterations ---------------------//
  for(int it = 1; it<=MAX_ITS; it++)
  { //printf("    Broyden: Iteration %d...\n", it);

    (obj->*func)(V);

    //if initial guess satisfies Accr, do not use broyden at all 
    if (it==1)
    {
      double MaxDiff = 0;
      for (int i=0; i<N; i++)
        if( abs( V[i] - Vinit[i] ) > MaxDiff ) 
           MaxDiff = abs( V[i] - Vinit[i] );
      if (MaxDiff < Accr) 
      {  printf("    !!! Initial guess good !!!\n");
         b = true;
         break;
      }
    }

    //save diff to compare it to new diff
    double OldDiff = B.CurrentDiff;
    
    int ErrCode = B.CalculateNew(V, it);
    if (ErrCode==1) { b=true; break; }
    else if (ErrCode==-1) { b=false; break; }
    
    
    //if Diff has increased by an order of magnitude that's a good sign 
    //Broyden is going further from the solution. Revert to initial guess
    if ( ( (B.CurrentDiff > 5.0 * OldDiff) and (it>3) ) 
         or 
         ( (B.CurrentDiff > 100.0) and (it>3) ) 
       )
    {  for (int i=0; i<N; i++) V[i] = Vinit[i];
       (obj->*func)(V);
       printf("    !-!!!! Broyden Stranded, Continuing with initial guess... !!!!!\n");
       break; 
    }
  }
  //--------------------------------------------//

  //turn off Broyden
  B.TurnOff();

  return b;
}


