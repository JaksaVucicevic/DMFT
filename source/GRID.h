#include <complex>

using namespace std;

namespace GridTypes
{
  const int LogLin = 0;
  const int Jaksa = 1;
  const int Linear = 2;
  const int MatsubaraLike = 3;
  const int LogLinEdge = 4;
}

class GRID
{
  private:
    //--grid parameters--//
    int GridType;
    int N;			//number of freq points
    int Nlog;			//number of freq points in logarithmic grid
    int Nlin;			//number of freq points in linear grid
    int Nedge;
    double omega_lin_max;	//freq cutoff
    double omega_max;		//maximum value of freq in the lograithmic grid
    double omega_min;		//lowest value of freq in the lograithmic grid
    double omega_edge;

    double domega_min;
    double domega_max;

    double* omega;

    void Defaults();
    
  public:
    GRID();
    GRID(int N, double omega_lin_max, bool OnlyPositive);
    GRID(int Nlog, int Nlin, double omega_lin_max, double omega_max, double omega_min);
    GRID(int Nlog, int Nlin, int Nedge, double omega_lin_max, double omega_max, double omega_min, double omega_edge);
    GRID(double domega_min, double domega_max, double omega_max, double omega_lin_max);
    GRID(const char* ParamsFN);
    ~GRID();
    
    int get_N() { return N; };
    double get_omega(int i, int Nlog, double omega_max,double omega_min);
    double get_omega_edge(int i);
    double get_omega(int i);
    double get_omega_lin_max() { return omega_lin_max; };
    double get_domega();
    double get_domega(double omega);
    void assign_omega(double* omega);
    
    //------routines--------//
    void KramarsKronig(complex<double> Y[]);
    complex<double> interpl(complex<double> X[], double om);
    double interpl(double X[], double om);
};
