#include <complex>
using namespace std;

class GRID;

class Result
{
  public:
    Result(GRID* grid);
    Result(const Result &result);
    ~Result();

    void Reset();
    void Reset(GRID* grid);
 
    GRID* grid;

    double n;
    double n0;
    double mu;
    double mu0;

    double* omega;		//omega grid
    double* fermi;		//fermi function
    double* Ap;			//spectral functions
    double* Am;
    double* P1;			//polarizations
    double* P2;
    complex<double>* SOCSigma;	//Second order contribution in sigma
    complex<double>* Sigma;	//Sigma interpolating between exact limiting cases
    complex<double>* G;		//Greens function on real axis
    complex<double>* Delta;	//Bath
    complex<double>* G0;	//auxillary Green's function
    double* DOS;		//quasi-particle density of states (typical DOS in TMT)
    double* NIDOS;		//non-interacting density of states
    double* DOSmed;		//medium DOS in TMT, can be used as an auxiallry DOS in other cases

    void PrintResult(const char* ResultFN);
    bool ReadFromFile(const char* ResultFN);
    void CopyFrom(const Result &result);
    void PrintModel(double U, double T);
    void PrintOnImagAxis(complex<double> * X, int M, double T, const char* FN);
    static void PrintCondOnImagAxis(int N, double* X, double* w, int M, double T, const char* FN);
    double get_ImGiw1(double T);
    
    double getIntegrand(double w, double T, double mu, double nu, double eps, complex<double> Sigma_nu,complex<double> Sigma_nu_plus_w);
    double Conductivity(double T, double mu, int Neps, int Nnu, const char * integrandFN = NULL, bool excludeSmallOmega = false);
    double TriangularConductivity(double T, int Nkx, int Nky, int Nnu, const char * integrandFN=NULL);
    double Conductivity(double w, double T, double mu, double Neps, double Nnu, const char * integrandFN = NULL);
    double NIConductivity(double T, double mu, int Neps, int Nnu, const char * integrandFN);
    void ChargeSusceptibility(double T, double &chi1, double &chi3);
    void PrintSpectralFunction(const char* FN, double (*eps)(double, double));
    void PrintFermiSurface(const char* FN, double (*eps)(double, double));

  private:
    void Initialize(GRID* grid);
    void ReleaseMemory();
};
