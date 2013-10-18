#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include "Input.h"
using namespace std;


Input::Input(const char* InputFN)
{
  SetInputFN(InputFN);
}

Input::~Input()
{
  
}

void Input::SetInputFN(const char* InputFN)
{ 
  this->InputFN.assign(InputFN);
  FILE* InputFile = fopen(this->InputFN.c_str(),"r");
  if ( InputFile==NULL )
    printf("-- WARNING -- Input: Input File does not exist!\n");
  else
  {  cout << "-- INFO -- Input: Input File name set to:" << this->InputFN << endl;
     fclose(InputFile);
  }
}
/*
int Input::ReadArray(int N, double* Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  int Error = 0;
  switch (N)
  {
    case 0: printf("-- ERROR -- Input: Arrays of 0 elements make no sense!\n"); return -1;
    case 1: Error = sscanf(line,"%le",                 &(Param[0]) ); break;
    case 2: Error = sscanf(line,"%le %le",             &(Param[0]), &(Param[1]) ); break;
    case 3: Error = sscanf(line,"%le %le %le",         &(Param[0]), &(Param[1]), &(Param[2]) ); break;
    case 4: Error = sscanf(line,"%le %le %le %le",     &(Param[0]), &(Param[1]), &(Param[2]), &(Param[3]) ); break;
    case 5: Error = sscanf(line,"%le %le %le %le %le", &(Param[0]), &(Param[1]), &(Param[2]), &(Param[3]), &(Param[4]) ); break;
    default: printf("-- ERROR -- Input: Arrays of more than 5 elements can not be read!\n"); return -1;
  }

  if ( Error == EOF ) 
  {  
     printf("-- ERROR -- Input: Param %s can not be read\n",ParamName);
     return -1;    
  }
  return 0;
}

int Input::ReadArray(int N, int* Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  int Error = 0;
  switch (N)
  {
    case 0: printf("-- ERROR -- Input: Arrays of 0 elements make no sense!\n"); return -1;
    case 1: Error = sscanf(line,"%d",             &(Param[0]) ); break;
    case 2: Error = sscanf(line,"%d %d",          &(Param[0]), &(Param[1]) ); break;
    case 3: Error = sscanf(line,"%d %d %d",       &(Param[0]), &(Param[1]), &(Param[2]) ); break;
    case 4: Error = sscanf(line,"%d %d %d %d",    &(Param[0]), &(Param[1]), &(Param[2]), &(Param[3]) ); break;
    case 5: Error = sscanf(line,"%d %d %d %d %d", &(Param[0]), &(Param[1]), &(Param[2]), &(Param[3]), &(Param[4]) ); break;
    default: printf("-- ERROR -- Input: Arrays of more than 5 elements can not be read!\n"); return -1;
  }

  if ( Error == EOF ) 
  {  
     printf("-- ERROR -- Input: Param %s can not be read\n",ParamName);
     return -1;    
  }
  return 0;
}*/

template <typename T> int Input::ReadArray(int N, T* Param,const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;
  stringstream ss;
  ss << line;
  for(int i = 0; i < N; i++)
    ss >> Param[i];
  return 0;
}

int Input::ReadArray(int N, double* Param, const char* ParamName)
{
  int Err;
  Err = ReadArray<double>(N, Param, ParamName);
  return Err;
}

int Input::ReadArray(int N, int* Param, const char* ParamName)
{
  int Err;
  Err = ReadArray<int>(N, Param, ParamName);
  return Err;
}

int Input::ReadParam(int& Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  if ( sscanf(line,"%d", &Param) == EOF ) 
  {  
     printf("-- ERROR -- Input: Param %s can not be read\n",ParamName);
     return -1;    
  }
  return 0;
}

int Input::ReadParam(double& Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  if ( sscanf(line,"%le", &Param) == EOF ) 
  {  
     printf("-- ERROR -- Input: Param %s can not be read\n",ParamName);
     return -1;    
  }
  return 0;
}

int Input::ReadParam(char& Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  Param = line[0]; 
  return 0;
}

int Input::ReadParam(char* Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  sscanf(line,"%s",Param);
  return 0;

}

int Input::ReadParam(bool& Param, const char* ParamName)
{
  char* line = ReadParam(ParamName);
  if (line==NULL) return -1;

  if (line[0]=='T')
    Param = true;
  else 
    if (line[0]=='F')
      Param = false;
    else
    {  printf("-- ERROR -- Input: Param %s can not be read\n",ParamName);
       return -1; 
    }
  return 0;
}


char* Input::ReadParam(const char* ParamName)
{ 
  FILE* InputFile = fopen(InputFN.c_str(),"r");
  if ( InputFile==NULL ) return NULL;
    
  char* line = new char[128];
  while ( fgets(line, 128, InputFile ) != NULL ) 
  { if (line[0]=='#') continue;
    string sline(line);
    if ( sline.find(ParamName) != string::npos )
      return line;
  }
  delete [] line;
  fclose (InputFile);
  printf("-- INFO -- Input: Param %s not found in Input File\n",ParamName);
  return NULL;
}
