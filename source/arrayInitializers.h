template<typename T> T**  Array2D(int N1, int N2)
{
  T** X = new T*[N1];
  for(int i=0; i<N1; i++)
  { X[i] = new T[N2];   
    for(int j=0; j<N2; j++)
      X[i][j] = 0.0;
  }
  return X;
}

template<typename T> T***  Array3D(int N1, int N2, int N3)
{
  T*** X = new T**[N1];
  for(int i=0; i<N1; i++)
    X[i] = Array2D<T>(N2,N3);
  return X;
}

template<typename T> T***  Array3D(int N1, int* N2, int n)
{
  T*** X = new T**[N1];
  for(int i=0; i<N1; i++)
    X[i] = Array2D<T>(N2[i]*n,N2[i]*n);   
  return X;
}

template<typename T> T****  Array4D(int N1, int N2, int* N3, int n)
{
  T **** X = new T***[N1];
  for(int i=0; i<N1; i++)
    X[i] = Array3D<T>(N2, N3[i]*n,N3[i]*n);   
  return X;
}

template<typename T> T***  Array4D(int N1, int N2, int N3, int N4)
{
  T **** X = new T***[N1];
  for(int i=0; i<N1; i++)
    X[i] = Array3D<T>( N2, N3, N4 );
  return X;
}

template<typename T> void  FreeArray2D(T** &X, int N1)
{
  for(int i=0; i<N1; i++)
    delete [] X[i];      
  delete [] X;
}

template<typename T> void  FreeArray3D(T*** &X, int N1, int N2)
{
  for(int i=0; i<N1; i++)
    FreeArray2D<T>(X[i], N2);
  delete [] X;
}

template<typename T> void  FreeArray3D(T*** &X, int N1, int* N2, int n)
{
  for(int i=0; i<N1; i++)
    FreeArray2D<T>(X[i], N2[i]*n);
  delete [] X;
}

template<typename T> void  FreeArray4D(T*** &X, int N1, int N2, int N3)
{
  for(int i=0; i<N1 i++)
    FreeArray3D<T>(X[i], N2, N3);
  delete [] X;
}

template<typename T> void  FreeArray4D(T**** &X, int N1, int N2, int* N3, int n)
{
  for(int i=0; i<N1; i++)
    FreeArray3D<T>(X[i], N2, N3[i]*n);
  delete [] X;
}

template<typename T> void Zeros2D(T** X, int N)
{
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)
    X[i][j] = (T) 0.0;    
}


template<typename T> void Zeros3D(T*** X, int N1, int* N2, int n)
{
  for(int i=0; i<N1; i++)
  for(int j=0; j<N2[i]*n; j++)
  for(int k=0; k<N2[i]*n; k++)
    X[i][j][k] = (T) 0.0;    
}

template<typename T> void Zeros3D(T*** X, int N1, int N2, int N3)
{
  for(int i=0; i<N1; i++)
  for(int j=0; j<N2; j++)
  for(int k=0; k<N3; k++)
    X[i][j][k] = (T) 0.0;    
}

template<typename T> void OneMinus(T** X, T** one_minus_X, int N)
{
  for(int i=0; i<N; i++)
  for(int j=0; j<N; j++)
    one_minus_X[i][j] = ((i==j)?(T)1.0:(T)0.0) - X[i][j];    
}



