#include <iostream>
#include "Loop.h"
using namespace std;

class SIAM;

class CHM: public Loop
{
  protected:

    void Defaults();
    string ParamsFN;

    double U;
    double T;
    double t;

    SIAM* siam;
    double SIAMeta;
    bool UseSmartSIAMeta;  

    bool UseBethe;
    int SiamNt;

    bool SIAMUseLatticeSpecificG;
    int LatticeType;
 
    virtual bool SolveSIAM();
    virtual void CalcDelta();  
   
    virtual void ReleaseMemory();

  public:
    CHM();
    CHM(const char* ParamsFN);   
    ~CHM();
  
    void SetParams(double U, double T, double t);
    void SetSIAMUseLatticeSpecificG(bool SIAMUseLatticeSpecificG);
    bool UseFixedMuSIAMRun;

    void SetUseBethe(bool UseBethe);
    void SetSIAMeta(double eta, bool UseSmartSIAMeta=false);

    double get_U() { return U; };
    double get_T() { return T; };
};


