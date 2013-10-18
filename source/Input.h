#include <iostream>
using namespace std;

class Input
{  
   public:
     Input(const char* InputFN);
     ~Input();

     void SetInputFN(const char* InputFN);   

     int ReadParam(double& Param, const char* ParamName);
     int ReadParam(int& Param, const char* ParamName);
     int ReadParam(char& Param, const char* ParamName);
     int ReadParam(char* Param, const char* ParamName);
     int ReadParam(bool& Param, const char* ParamName); 
     int ReadArray(int N, double* Param, const char* ParamName);
     int ReadArray(int N, int* Param, const char* ParamName);

   private:
     string InputFN;
     template <typename T> int ReadArray(int N, T* Param,const char* ParamName);
     char* ReadParam(const char* ParamName);
};
