#include <complex>

using namespace std;

namespace IBZtypes
{
  const int SquareLattice = 0;
  const int TriangularLattice = 1;
}

class IBZ
{
  private:
    int IBZtype;
    
    int Nx;
    int Ny;
    int Ntotal;

    double* kx;
    double* ky;

    double a; //lattice spacing;
    //double b;

    double get_kx(int i);
    double get_ky(int j);
    double get_epsilon(double kx, double ky);
    double get_velocity(double kx, double ky);
    void Initialize();
    void ReleaseMemory();

    double t;
  public:
    IBZ(int IBZtype, int Nx, int Ny, double a=1.0);
    ~IBZ();

    void SetLatticeSpacing(double a);  

    double** epsilon;
    double** velocity;
    complex<double>** summand;
  
    complex<double> sum();
    void PrintToFile(const char* FN);
};
