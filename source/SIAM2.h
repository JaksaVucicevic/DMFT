//**********************************************************//
//              SIAM at Arbitrary Filling                   //
//                    by Jaksha Vuchichevicc                //
//**********************************************************//

#include <iostream>
#include <complex>

class Result;
class GRID;

using namespace std;

namespace FormulasForG
{
  const int LatticeSpecific = 0;
  const int NoLattice = 1;
  const int Integral = 2;

}

//======================= SIAM Class ==========================================//

class SIAM2
{
  private:

    void Defaults();
    string ParamsFN;

    Result* r;

    //--impurity parameters--//
    double U;			//on-site repulsion
    double T;			//temperature
    double epsilon;		//impurity energy level

    //---BROADENING---//
    double eta;   

    //----lattice---------//
   
    //-- Broyden solver options--//
    double Accr;
    int MAX_ITS;  

    //--MPT Higher order correlations--//
    double MPT_B;
    double MPT_B0;
    
    //--storage arrays--//
    GRID* grid;
    int N;

    //--get functions--//
    double get_fermi(int i);
    double get_n(complex<double> X[]);

    //--get procedures--//
    void get_fermi();
    void get_G0();
    void get_As();
    void get_Ps();  
    void get_SOCSigma();
    double get_MPT_B();
    double get_MPT_B0();
    double get_b();
    double get_a();
    void get_Sigma();
    void get_G();

    bool ClipOff(complex<double> &X);
    bool Clipped;

    //--imaginary axis--// 
    double MatsFreq(int n);

  public:
    //------ OPTIONS -------//
    bool UseMPT_Bs;		//if true program uses MPT higher coerrelations B and B0
    bool CheckSpectralWeight;   //if true program prints out spectral weights of G and G0 after each iteration
    void SetBroydenParameters(int MAX_ITS, double Accr);
    void SetBroadening(double eta);
    void SetDOStype_CHM(int DOStype, double t, const char* FileName ="");

    void SetT(double T);
    void SetU(double U);
    void SetEpsilon(double epsilon);
    void SetUTepsilon(double U, double T, double epsilon);

    int FORMULA_FOR_G;
    int LatticeType;
    double t;
    void SetUseLatticeSpecificG(double t, int LatticeType);
    bool HalfFilling;

    //--Constructors/destructors--//
    SIAM2();  
    SIAM2(const char* ParamsFN);
    ~SIAM2();
    
    //get G on inamginary axis
    void GetGfOnImagAxis(int Nmax, complex<double>* G_out);
    
    //--------RUN SIAM--------//
    
    bool Run(Result* r); 

    //--print out routines--//
    void PrintModel();

  //---- FRIENDS -----//
  //function that will be calling private member functions in case of solving (systems of) equations
  friend bool UseBroyden(int, int, double, void (SIAM2::*)(complex<double>*), SIAM2*, complex<double>*);
   
};
