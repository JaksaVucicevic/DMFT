#include "SIAM.h"
#include "Mixer.h"

using namespace std;

/*************************************/
//       Clean Hubbard Model         //
/*************************************/


class CHM
{
  private:
    bool Initialized;
    
    GRID * grid;
    int N;
    double * omega;
    complex<double>* LastDelta;

    double Accr;

    int MAX_ITS;
    int NtoMix;
    const int * Coefs;

    bool PrintIntermediate;
    bool HaltOnIterations;
   
    bool StartFromInsulator;
    bool StartFromPrevious;
    bool UseBethe;

    bool UseBroyden;
    bool ForceBroyden;
    double BroydenStartDiff;
    
    double SIAMeta;
  
  public:
    CHM();
    ~CHM();

    void SetGrid(GRID * grid);
    void SetMaxItsAccrAndCoefs(int MAX_ITS, double Accr, int NtoMix, const int * Coefs);
    void SetPrintIntermediate(bool PrintINtermediate);
    void SetHaltOnIterations(bool HaltOnIterations);
    void SetBroyden(bool UseBroyden, bool ForceBroyden, double BroydenStartDiff);
    void SetStartFromInsulator(bool StartFromInsulator);
    void SetStartFromPrevious(bool StartFromPrevious);
    void SetUseBethe(bool UseBethe);
    void SetSIAMeta(double eta);

    void Run(double n, double U, double T, double t, double &mu, int DOStype, const char* FNDos = "");
    int IterationsMade;
    void InitDeltaFromFile(const char* FNDelta, int Nlog, int Nlin, double omega_lin_max, double omega_max, double omega_min); 
};


