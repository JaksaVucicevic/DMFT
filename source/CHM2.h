#include <iostream>
#include "Loop.h"
using namespace std;

class SIAM2;

class CHM2: public Loop
{
  protected:

    void Defaults();
    string ParamsFN;

    double U;
    double T;
    double t;

    
    double SIAMeta;
    bool UseBethe;
    int SiamNt;

    bool SIAMUseLatticeSpecificG;
    int LatticeType;
 
    virtual bool SolveSIAM();
    virtual void CalcDelta();  
   
    virtual void ReleaseMemory();

  public:
    CHM2();
    CHM2(const char* ParamsFN);   
    ~CHM2();
     
    SIAM2* siam;
    void SetParams(double U, double T, double t);
    void SetSIAMUseLatticeSpecificG(bool SIAMUseLatticeSpecificG);

    void SetSIAMeta(double eta);

    double get_U() { return U; };
    double get_T() { return T; };
};


