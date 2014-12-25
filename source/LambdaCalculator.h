#include <complex>

using namespace std;

#define MAX_LAMBDAS 10

class LambdaCalculator
{
  private:
    void Defaults();
    void Initialize();
    void ReleaseMemory();
    bool Initialized;
    //------iteration counter-------//
    int counter;

    //--------- input X ---------//
    int N;
    double* rX[3];			//input X
    complex<double>* cXw[3];		//input X
    
    int offset;					//used only for the input X. not used with continued X				
    int Nlambdas;				//number of different lambdas to be calculated (max 10)
    int Nfreqs[MAX_LAMBDAS];			//the number of freqs used for calculation of eachj lambda
    int Ndiffs;					//number of diffs to be calculated (max 10)
    int is[MAX_LAMBDAS];					//the index of X to be used in diff calculation. if -l, maximum difference is found

    //--------- continued X ---------//
    int M;				//number of matsubara freq to be saved
    complex<double>* cXiw[3];		//continued X
    double* omega; 			//used only in continuation to imag axis.
    double* iw;
    double T;  				//used only in continuation to imag axis.

    int continued_Nlambdas;				//number of different lambdas to be calculated (max 10)
    int continued_Nfreqs[MAX_LAMBDAS];				//the number of freqs used for calculation of eachj lambda
    int continued_Ndiffs;					//number of diffs to be calculated (max 10)
    int continued_is[MAX_LAMBDAS];					//the index of X to be used in diff calculation. if -l, maximum difference is found

    //--------- methods --------------//
    void PrintOut(); 
    void ShiftAddAndContinue(complex<double>* cX);

    //----------options----------------//
    bool DoContinued;				// calculate and output  lambdas and diffs from continued X
    bool DoOutput;				// if set to false, only returned lambda is calculated and no output files are created
    bool DoPrintOutContinuedX;
    bool SmartDiff;				//if set to false only smart diff is calculated

    char OutputFN[300];				// output options

  public:
    LambdaCalculator();
    LambdaCalculator(int N, int M, double T);
    LambdaCalculator(const char* paramsFN);
    ~LambdaCalculator();

    //options
    void SetN(int N);
    void SetT(double T);
    void SetNMT(int N, int M, double T);
    void SetOmega(double* omega);
    void SetOffset(int offset);
    void SetLambdasPrinted(int Nlambdas, const int* Nfreqs);
    void SetContinuedLambdasPrinted(int continued_Nlambdas, const int*continued_Nfreqs);
    void SetDiffsPrinted(int Ndiffs, const int* is);
    void SetContinuedDiffsPrinted(int continued_Ndiffs, const int* continued_is);
    void SetOutputFileName(const char* FN);	//only the parameters extnesion (e.g. ".U2.500.T0.050.W1.000.n0.500"). the part "diffs", "lambdas" or "lambdas.continued." is automaticaly inserted at the beginning
    void SetDoOutput(bool DoOutput);
    void ResetCounter(int counter=0);
     
    //====================================//
    double CalculateLambda(double* rX);		// continued X not applicable when this is called
    double CalculateLambda(complex<double>* cX);

    double lambdas[MAX_LAMBDAS];
    double diffs[MAX_LAMBDAS];
    double smart_diff;
    double continued_lambdas[MAX_LAMBDAS];
    double continued_diffs[MAX_LAMBDAS];
    double continued_smart_diff;
    //====================================//

    double lambda_history[1000];  
    double best_lambda;

    //actual code    
    static void ContinueToImagAxis(int N, complex<double>* Xw, double* w, int M, complex<double>* Xiw, double* iw, double T=-1.0, const char* outputFN=NULL);
    static double CalcLambda(int offset, int Nfreq, int N, complex<double>* X0, complex<double>* X1, complex<double>* X2);
                        //use Nfreq frequences starting form X[offset]. if  Nfreq==-1, use all N-offset frequences. 
    static double CalcDiff(int offset, int i, int N, complex<double>* X1, complex<double>* X2 );	//which freq X[i] is used for calculating diff
    static double CalcSmartDiff(int offset, int N, complex<double>* X1, complex<double>* X2 );
    static double FindBestLambda(int counter, int n, const double* lambda_history);
};
