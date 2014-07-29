//**********************************************************//
//              SIAM at Arbitrary Filling                   //
//                    by Jaksha Vuchichevicc                //
//**********************************************************//

#include <iostream>
#include <complex>

class Result;
class GRID;

using namespace std;

//======================= SIAM Class ==========================================//

class SIAM
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
    
    //--bath parameters--// 
    int DOStype_CHM;		//used in RunCHM for calculation of G
    double t_CHM;
    double mu;			//global chemical potential
    double mu0;			//fictious chemical potential
    bool isBethe;		//Set this to true when using bethe lattice specific self-consistency to use simplified expression for G
    bool GfromDelta;            // affects how G is calculated within SolveSiam function which is called within Run. If true, G = 1/( w+mu-Sigma-Delta )

    //----lattice---------//
    bool UseLatticeSpecificG;
    int LatticeType;
    double t;
    
    //--don't touch this---//
    bool SymmetricCase;
    bool HalfFilling;

    //-- Broyden solver options--//
    double Accr;
    int MAX_ITS; 

    //-- mu0 search --//
    bool UseBroydenFormu0;	//set this to false if only Amoeba method is tu be used
    int max_tries;		//number of tries (with different initial guesses) of broyden search before Amoeba is used
                                //set to a large number (say 100) if MPT corrections are used. Amoeba does only the mu0 search
    				

    //--MPT Higher order corrections--//
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
    void get_G0(complex<double>* V);
    void get_As();
    void get_Ps();  
    void get_SOCSigma();
    double get_MPT_B();
    double get_MPT_B0();
    double get_b();
    void get_Sigma();
    void get_G();
    void get_G(complex<double>* V); //used by broyden in solving systems of equations
    void get_G_CHM();
    void get_G_CHM(complex<double>* V); //used by broyden in solving systems of equations

    bool ClipOff(complex<double> &X);
    bool Clipped;

    //--imaginary axis--// 
    double MatsFreq(int n);

    //--- SIAM solver ---//
    void SolveSiam(complex<double>* V);
    void Amoeba(double accr, complex<double>* V);	//amoeba method for mu0 search. not applicable when MPT corrections are used (TODO: generalize this method)
    void Amoeba_CHM(double accr, complex<double>* V);

    double AmoebaScanStart;	//before amoeba starts, the equation is solved roughly (with accuracy AmobeScanStep) by scanning from AmoebaScanStart to AmoebaScanEnd.
    double AmoebaScanEnd; 	//make sure AmoebaScanStart and AmoebaScanEnd are far enough apart (when U or W is large).
    double AmoebaScanStep;
    int AmoebaMaxIts;		//maximum number of Amoeba iterations
    bool AmoebaForceScanAndPrintOut;	//output n, n0, n-n0 and result when scanning for mu0 candidate
  
  public:
    //------ OPTIONS -------//
    bool UseMPT_Bs;		//if true program uses MPT higher coerrelations B and B0
    bool CheckSpectralWeight;   //if true program prints out spectral weights of G and G0 after each iteration
    void SetBroydenParameters(int MAX_ITS, double Accr);
    void SetBroadening(double eta);
    void SetDOStype_CHM(int DOStype, double t, const char* FileName ="");
    void SetIsBethe(bool isBethe);
    void SetT(double T);
    void SetU(double U);
    void SetEpsilon(double epsilon);
    void SetUTepsilon(double U, double T, double epsilon);
    void SetAmoebaParams(double AmoebaScanStart, double AmoebaScanEnd, double AmoebaScanStep);

    void SetUseLatticeSpecificG(bool UseLatticeSpecificG, double t, int LatticeType);

    //--Constructors/destructors--//
    SIAM();  
    SIAM(const char* ParamsFN);
    ~SIAM();
    
    //get G on inamginary axis
    void GetGfOnImagAxis(int Nmax, complex<double>* G_out);

    static void get_G_from_Sigma(Result* r);
    
    //--------RUN SIAM--------//
    
    bool Run(Result* r); 
    bool Run_CHM(Result* r); 

    //--print out routines--//
    void PrintModel();

  //---- FRIENDS -----//
  //function that will be calling private member functions in case of solving (systems of) equations
  friend bool UseBroyden(int, int, double, void (SIAM::*)(complex<double>*), SIAM*, complex<double>*);
   
};
