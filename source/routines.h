#include <complex>
#include <vector>

using namespace std;

//======================== CONSTANTS ===============================//

const double pi = 3.14159265358979323846;
const double e = 2.71;
const complex<double> ii = complex<double>(0.0,1.0);

//======================= ROUTINES ==================================//

double sign(double x);
int int_sign(double x);

double sqr(double x);
int pow(int base, int exp);
complex<double> sqr(complex<double> x);
/*double abs(double x);*/
/*double abs(complex<double> x);*/

//--- integral ---//

double TrapezIntegral(int N, double Y[], double X[]);
complex<double> TrapezIntegral(int N, complex<double> Y[], double X[]);
double TrapezIntegralMP(int N, double Y[], double X[]);
complex<double> TrapezIntegralMP(int N, complex<double> Y[], double X[]);
double SmartIntegral(int N, double Y[], double X[], double x0, double Dx, double (*f)(double), int Nextra, double x_min, const char* FN = NULL );
//double TrapezIntegral(std::vector< double > Y, std::vector<double> X);
complex<double> TrapezIntegral(std::vector< complex<double> > Y, std::vector<double> X);
double EllipticIntegralFirstKind(double x);
double SI(double x);
complex<double> EllipticIntegralFirstKind(complex<double> x);
double interpl(int N, double* Y, double* X, double x);

//======================== IO =======================================//

void PrintFunc(const char* FileName, int N, int M, double** Y, double* X);
void PrintFunc(const char* FileName, int N, complex<double>* Y, double* X);
void PrintFunc(const char* FileName, int N, complex<double>* Y);
void PrintFunc(const char* FileName, int N, double* Y);
void PrintFunc(const char* FileName, int N, double* Y, double* X);
void PrintFunc(const char* FileName, std::vector< complex<double> > Y, std::vector<double> X);
void PrintFunc3D(const char* FileName, int N, complex<double>** Y, double* X);
void ReadFunc(const char* FileName, int &N, int &M, double** &X);
//void ReadFunc(const char* FileName, int &N, complex<double>* &Y, double* &X, bool PurelyReal=false);
//void ReadFunc(const char* FileName, int &N, double* &Y, double* &X);
void ReadFunc(const char* FileName, int N, double* Y, double* X);
bool FileExists(const char* FN);
//===================vectors and matrices=============================//

void MultiplyByMatrix(int N, double* v, double** m);
void CreateRotationMatrix(double** m, int N, double angle, int* plane);
void RotateVector(int N, double* v, double angle, int* plane);
void RotateVector2D(double* v, double angle);
void InvertMatrix(int N, double** A, double** invA, double &det);

//==================== DOSes and Init Deltas ========================//

namespace DOStypes
{
  const int SemiCircle = 0;
  const int Gaussian = 1;
  const int Insulator = 2;
  const int SquareLattice = 3;
  const int CubicLattice = 4;
  const int FromFileSymmetric = 5;
  const int FromFile = 6;
  const int Uniform = 7;
  //add more if needed
}

double DOS(int DOStype, double t, double om, double U=0.0);

void WriteCubicDosToFile();
void ReadDosFromFile(const char* FN, int N, double* omega, double* DOS);

void get_G_from_DOS(int DOStype, double t, int N, double* omega, complex<double>* G, double eta, double U=0.0);
void get_Sigma_from_G(double t, int N, double* omega, complex<double>* G, complex<double>* Sigma);

complex<double> LS_get_G(int DOStype, double t, complex<double> com);
void InitG(int DOStype, double t, int N, double* omega, complex<double>* G);
void InitG(int DOStype, double t, int N, complex<double>* omega, complex<double>* G);

void InitDOS(int DOStype, double t, int N, double* omega, double* dos, double U=0.0);
void InitDelta(int DOStype, 
               int N, 
               double V, 
               double mu, 
               double eta, 
               double t,
               double* omega,  
               complex<double>* Delta,
               const char* FNDos = "");

void InitDeltaFromSIAMOutputFile(const char* FN, int N, double* omega, complex<double>* Delta);

void InitInsulatorDelta(double U,
                        int N,
                        double t,
                        double* omega,
                        complex<double>* Delta);
                                                                                                
