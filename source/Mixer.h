#include <cstdlib>
#include "routines.h" 

template <class T> //T may be  float, double, complex<double>
class Mixer
{
  private:
    bool Initialized;
    
    int N; //length of solution array
    int M; //number of solution to save
    int Counter; //counts iterations    
    const int * Coefs; //mixing coefficients
    T** X; //saved solutions
    double Accr; //required accuracy
    
    bool CheckConvergence(); //returns true if solutions converged
    void AddSolution(T* Solution);
    void MixSolutions();
    void ShiftSolutions();
    void OutputSolution(T * Solution);
    
    void ReleaseMemory();
  public:
    Mixer(int N, int M, const int * Coefs, double Accr);
    Mixer();
    ~Mixer();
    
    //------ Call this to mix solutions and check convergence ----//
    bool Mix(T* Solution); //returns true if solutions converged
   
    double CurrentDiff;
     //initializes Mixer for M N-long Solutions that will be mixed with Coefs until Accr reached
    void Initialize(int N, int M, const int * Coefs, double Accr);
    void Reset();
};

template <class T>
Mixer<T>::Mixer(int N, int M, const int * Coefs, double Accr)
{
  Initialized = false;
  Initialize(N, M, Coefs, Accr);
}

template <class T>
Mixer<T>::Mixer()
{
  Initialized = false;
}

template <class T>
Mixer<T>::~Mixer()
{
  ReleaseMemory();
}

template <class T>
void Mixer<T>::Initialize(int N, int M, const int * Coefs, double Accr)
{
  if (Initialized) ReleaseMemory();
  this->N = N;
  this->M = M;
  X = new T*[M];
  for (int n=0; n<M; n++) X[n] = new T[N];
  this->Coefs =  Coefs;
  this->Accr = Accr;
  Counter = 0;
  CurrentDiff = 1e-2;
  Initialized = true;
}

template <class T>
void Mixer<T>::ReleaseMemory()
{
  for (int n=0; n<M; n++) delete [] X[n];
  delete [] X;
}

template <class T>
void Mixer<T>::Reset()
{
  Initialize(N,M,Coefs,Accr);
}

template <class T>
void Mixer<T>::AddSolution(T* Solution)
{
  for(int i=0; i<N; i++) X[0][i] = Solution[i];
}

template <class T>
void Mixer<T>::ShiftSolutions()
{
  for(int i=0; i<N; i++)
    for(int n=M-1; n>0; n--) 
      X[n][i]=X[n-1][i]; //shift old solutions
}

template <class T>
void Mixer<T>::MixSolutions()
{
  T sum;
  int nds = (Counter<M) ? Counter-1 : M-1;
  int denom = 0;
  for(int n=0; n<=nds ; n++)
    denom += Coefs[n];
  
  for(int i=0; i<N ; i++)
  { 
    sum = (T) 0.0;  
    for(int n=0; n<=nds ; n++)
      sum += (T)Coefs[n] * X[n][i];

    X[0][i] = sum / (T)(denom);
  }
}

template <class T>
void Mixer<T>::OutputSolution(T * Solution)
{
  for(int i=0; i<N; i++) Solution[i] = X[0][i];
}

template <class T>
bool Mixer<T>::Mix(T* Solution)
{
  bool b;
  if (!Initialized) { printf("----Mixer: ERROR: Mixer Not Initialized !!!! Exiting to system... \n"); exit(1); }
  Counter++;
  if (Counter>1) 
  { ShiftSolutions();
    AddSolution(Solution);
    b = CheckConvergence();
    if (!b)
    {  MixSolutions();
       OutputSolution(Solution);
    }
  }
  else AddSolution(Solution);
  return b;
}

template <class T>
bool Mixer<T>::CheckConvergence()
{ 
  double MaxDiff=0; 
  for(int i=0; i<N; i++)
    if ( abs( X[0][i] - X[1][i] ) > MaxDiff )
      MaxDiff =  abs( X[0][i] - X[1][i] ) ;
  CurrentDiff = MaxDiff;
  printf("--- Mixer: Diff[%d] = %le ---\n", Counter, MaxDiff);

/*
  FILE* diffsFile = fopen("diffs","a");
  fprintf(diffsFile,"%le\n", MaxDiff);
  fclose(diffsFile);
*/
  if (MaxDiff < Accr) printf("--- Mixer: CONVERGED !!!\n");
  return  (MaxDiff < Accr);
}

