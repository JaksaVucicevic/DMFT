#include <cstdio>
#include "GRID.h"
#include "routines.h"
#include "Input.h"
#include <omp.h>

//==================== Constructors destruictors ==================//

void GRID::Defaults()
{
  GridType = GridTypes::LogLin;
  Nlog = 1500; //default 800
  Nlin = 1500; //default 570
  N = Nlog+Nlin;
  omega_lin_max = 6.0; //default 8.0
  omega_max = 0.3;   //default 1.0
  omega_min = 1e-10; //default 1e-6
}

GRID::GRID()
{
  Defaults();
}

GRID::GRID(int N, double omega_lin_max, bool OnlyPositive)
{
  if (!OnlyPositive)
    GridType = GridTypes::Linear;
  else 
    GridType = GridTypes::MatsubaraLike;
  this->N = N;
  this->omega_lin_max = omega_lin_max;

}

GRID::GRID(int Nlog, int Nlin, double omega_lin_max, double omega_max, double omega_min)
{
  GridType = GridTypes::LogLin;
  this->Nlog = Nlog;
  this->Nlin = Nlin;
  this->N = Nlog+Nlin;
  this->omega_lin_max = omega_lin_max;
  this->omega_max = omega_max;
  this->omega_min = omega_min;
}

GRID::GRID(int Nlog, int Nlin, int Nedge, double omega_lin_max, double omega_max, double omega_min, double omega_edge)
{
  GridType = GridTypes::LogLinEdge;
  this->Nlog = Nlog;
  this->Nlin = Nlin;
  this->Nedge = Nedge;
  this->N = Nlog+Nlin+Nedge;
  this->omega_lin_max = omega_lin_max;
  this->omega_max = omega_max;
  this->omega_min = omega_min;
  this->omega_edge = omega_edge;
}


GRID::GRID(double domega_min, double domega_max, double omega_max, double omega_lin_max)
{
  GridType = GridTypes::Jaksa;
  this->omega_min = domega_min;
  this->omega_lin_max = omega_lin_max;
  this->omega_max = omega_max;
  this->domega_max = domega_max;
  this->domega_min = domega_min;
}



GRID::GRID(const char* ParamsFN)
{
  Defaults();

  Input input(ParamsFN);
  input.ReadParam(GridType,"GRID::GridType");
  switch (GridType)
  { case GridTypes::LogLin :
    {
      if (input.ReadParam(Nlog,"GRID::Nlog")==-1) printf("-- INFO -- GRID: Continuing with default value...\n");
      if (input.ReadParam(Nlin,"GRID::Nlin")==-1) printf("-- INFO -- GRID: Continuing with default value...\n");
      N = Nlog + Nlin;
      if (input.ReadParam(omega_lin_max,"GRID::omega_lin_max")==-1) printf("-- INFO -- GRID: Continuing with default value...\n");
      if (input.ReadParam(omega_max,"GRID::omega_max")==-1) printf("-- INFO -- GRID: Continuing with default value...\n");
      if (input.ReadParam(omega_min,"GRID::omega_min")==-1) printf("-- INFO -- GRID: Continuing with default value...\n");
    } break;
    case GridTypes::Jaksa :
    {
      // TODO
    } break;
    case GridTypes::Linear :
    case GridTypes::MatsubaraLike :
    { if (input.ReadParam(N,"GRID::N")==-1) 
        printf("-- INFO -- GRID: Continuing with default value...\n");
      else
        printf("-- INFO -- GRID: N set to %d\n",N);
      if (input.ReadParam(omega_lin_max,"GRID::omega_lin_max")==-1) 
        printf("-- INFO -- GRID: Continuing with default value...\n");
      else
        printf("-- INFO -- GRID: omega_lin_max set to %f\n",omega_lin_max);
 
      // TODO
    } break;
  }
}

GRID::~GRID()
{

}
//======================= Initializers =============================//

double GRID::get_omega(int i, int Nlog, double omega_max,double omega_min)
{ 
  int sgn1 = (i >= Nlog/2) ? 1 : -1;
  i = (i >= Nlog/2) ? i-Nlog/2 : (Nlog/2-1)-i;
  return   sgn1 * ( exp ( log(omega_min)
                          + (double)i / (Nlog/2.0-1.0)
                            * log(omega_max/omega_min) ) ); 
}

double GRID::get_omega_edge(int i)
{
  if ((i<Nedge/2)or(i>=Nedge/2 + Nlog + Nlin))
  {
    if (i<Nedge/2)
      return get_omega(i+Nedge/2, Nedge, omega_edge-omega_lin_max-(get_omega(1)-get_omega(0)), omega_min) - omega_edge;
    else
      return get_omega(i-(N-Nedge/2), Nedge, omega_edge-omega_lin_max-(get_omega(1)-get_omega(0)), omega_min) + omega_edge;
  }
  else
    return get_omega(i-Nedge/2);

}

double GRID::get_omega(int i)
{
  if ((i<Nlin/2)||(i>= Nlin/2 + Nlog))
  {
    return (i<Nlin/2) ? -omega_lin_max + i*(omega_lin_max-omega_max)/(Nlin/2)
                      : omega_max + (i - Nlin/2 - Nlog + 1)*(omega_lin_max-omega_max)/(Nlin/2);
  }
  else
  { i-=Nlin/2;
    int sgn1 = (i >= Nlog/2) ? 1 : -1;
    i = (i >= Nlog/2) ? i-Nlog/2 : (Nlog/2-1)-i;
    return   sgn1 * ( exp ( log(omega_min)
                            + (double)i / (Nlog/2.0-1.0)
                              * log(omega_max/omega_min) ) );
  }
}


double GRID::get_domega()
{
  switch (GridType)
  { case GridTypes::MatsubaraLike:
      return 2.0 * omega_lin_max / ( 2.0 * (N-1.0) + 1.0 ) ;
      break;
    default: return 0;
  }
}


double GRID::get_domega(double omega)
{
  if (omega>=omega_max) return domega_max;
  else return (domega_max-domega_min)/omega_max * omega + domega_min;
}

void GRID::assign_omega(double* omega)
{ 
  this->omega = omega;
  switch (GridType)
  {
    case GridTypes::LogLin:
    {
      for (int i=0; i<N; i++) omega[i] = get_omega(i);
    } break;
    case GridTypes::Jaksa:
    {   
      int count=0;
      for(double w=domega_min/2.0; w<omega_lin_max; w+=get_domega(w))
        count++;
      N=2*count;
      for(int i=N/2; i<N; i++)
      {  omega[i] = (i==N/2) ? domega_min/2.0 : (omega[i-1] + get_domega(omega[i-1]));
         omega[N-1-i] = - omega[i];
      }
      omega_lin_max = omega[N-1];
      printf(">>>>>> GRID: N=%d, omega_lin_max=%.6f\n", N, omega_lin_max);
    } break;
    case GridTypes::Linear:
    { double domega = 2.0 * omega_lin_max / (N-1);
      for (int i=0; i<N; i++) omega[i] = - omega_lin_max + domega*i;
    } break;
    case GridTypes::MatsubaraLike:
    { double domega = get_domega();
      for (int i=0; i<N; i++) omega[i] = domega * (i + 0.5) ;
    } break;
    case GridTypes::LogLinEdge:
    {
      for (int i=0; i<N; i++) omega[i] = get_omega_edge(i);
    } break;
    
  }
}

//================================ routines ===============================//

void GRID::KramarsKronig(complex<double> Y[])
{
 
  if (omega == NULL) 
  {
    printf("-- Error -- GRID: KramarsKronig: No omega array assigned");
    return;
  }

  //int k;

  double** s = new double*[N];  
  #pragma omp parallel for
  for (int i=0; i<N; i++)
  { 

    double y = imag(Y[i]);
    double LogTerm = ( (i==0) || (i==N-1) ) 
                    ? 0.0
                    : y * log( (omega_lin_max-omega[i])
                              /(omega[i]+omega_lin_max) );

    s[i] = new double[N];
    for (int j=0; j<N; j++)
    { 
      if (i==j)
        s[i][j] =  imag(   Y[ j + ( (j < N-1) ? 1 : 0 ) ]
                      - Y[ j - ( (j > 0)   ? 1 : 0 ) ] )   
               / (  omega[ j + ( (j < N-1) ? 1 : 0 ) ]
                  - omega[ j - ( (j > 0)   ? 1 : 0 ) ] );
      else
        s[i][j] = ( imag(Y[j]) - y) 
               / ( omega[i]-omega[j] );
    }                     
    Y[i] = complex<double>( - ( TrapezIntegral(N, s[i], omega) - LogTerm )/pi , y);

    delete [] s[i];
  }
  delete [] s;
}
/*
void GRID::KramarsKronig(complex<double> Y[])
{
  if (omega == NULL) 
  {
    printf("-- Error -- GRID: KramarsKronig: No omega array assigned");
    return;
  }

  int k;
  for (int i=0; i<N; i++)
  { 
    double* s = new double[N];  
    double y = imag(Y[i]);
    double LogTerm = ( (i==0) || (i==N-1) ) 
                    ? 0.0
                    : y * log( (omega_lin_max-omega[i])
                              /(omega[i]+omega_lin_max) );

    for (int j=0; j<N; j++)
    { 
      if (i==j)
        s[j] =  imag(   Y[ j + ( (j < N-1) ? 1 : 0 ) ]
                      - Y[ j - ( (j > 0)   ? 1 : 0 ) ] )   
               / (  omega[ j + ( (j < N-1) ? 1 : 0 ) ]
                  - omega[ j - ( (j > 0)   ? 1 : 0 ) ] );
      else
        s[j] = ( imag(Y[j]) - y) 
               / ( omega[i]-omega[j] );
    }                     
    Y[i] = complex<double>( - ( TrapezIntegral(N, s, omega) - LogTerm )/pi , y);

    delete [] s;
  }
}
*/

double GRID::interpl(double X[], double om)
{   
  if (omega == NULL) 
  {
    printf("-- Error -- GRID: interpl: No omega array assigned");
    return 0;
  }
  
  switch (GridType)
  {  //---------- LOG LIN --------------//
     case GridTypes::LogLin:
     {
        if (abs(om) > omega_lin_max) 
          return 0.0;
        else
        {

          if ( (abs(om) > omega_max) && (abs(om) <= omega_lin_max) )  
          {  if (om > 0.0)
            {
              double dc = (om - omega_max)/(omega_lin_max-omega_max)*Nlin/2;
              int c = (int) dc;
              int k = c + Nlin/2 + Nlog - 1;
              double m = 1.0;
              if ((Nlog==0)and(k==N/2-1)) m=0.5;
              return X[k]+(X[k+1]-X[k])*(dc - (double)c)*m + (1.0-m)*(X[k+1]-X[k]);
            }
            else
            {
              double dc = (om + omega_lin_max)/(omega_lin_max-omega_max)*Nlin/2;
              int c = (int) dc;
              int k = c;
              double m = 1.0;
              if ((Nlog==0)and(k==N/2-1)) m=0.5;
              return X[k]+(X[k+1]-X[k])*(dc - (double)c)*m;
            }
          }
   
          if (abs(om) <= omega_min)
          {  if (omega_min==0) 
               return  0.5*( X[N/2-1] + X[N/2] );
             else 
               return X[N/2-1] + (om-(-omega_min))/(2*omega_min)*(X[N/2]-X[N/2-1]);
          }

          if ( (abs(om) > omega_min) and (abs(om) <= omega_max))
          {
            //if (Nlog==0) return 0.5*( X[N/2-1] + X[N/2] );
            double dc = (Nlog/2.0-1.0) * log(abs(om)/omega_min)
                                  / log(omega_max/omega_min);
            int c = (om>0) ? (int) dc + Nlog/2 : Nlog/2 - 2 - (int) dc;         
            int k = c + Nlin/2;
            double t = (om - omega[k])/(omega[k+1]-omega[k]);
    
            return  X[k]+(X[k+1]-X[k])*t;
          }
        }
     } break;
     //-------------- Jaksa Type ---------------//
     case GridTypes::Jaksa:
     {
        if (abs(om) > omega_lin_max) 
          return 0.0;
        else
        { int i;
          for (i = (om<0) ? 0 : N/2; i<N; i++)
            if (omega[i]>om) break;
          return X[i-1] + (X[i]-X[i-1]) / (omega[i]-omega[i-1]) 
                          * (om-omega[i-1]);
        }    

     } break;
     //------------- PURE LINEAR ----------------//
     case GridTypes::Linear:
     { if (abs(om) > omega_lin_max) 
          return 0.0;
       double domega = omega_lin_max / (N-1);
       int i = (int) (om+omega_lin_max)/domega ; 
       double c = om - omega[i];
       return c * (X[i+1]-X[i]) / domega;
     } break;
  }
  return 0;
}

complex<double> GRID::interpl(complex<double> X[], double om)
{
  double* Re = new double[N];
  double* Im = new double[N];
  for (int i=0; i<N; i++)
  {
    Re[i] = real(X[i]);
    Im[i] = imag(X[i]);
  }
  delete [] Re;
  delete [] Im;
  return complex<double>(interpl(Re,om),interpl(Im,om));
}
