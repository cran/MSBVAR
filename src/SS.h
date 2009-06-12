#ifndef SS_H
#define SS_H
#define UINT uint32_t
#define NBITS 30
#define NBITSF 30.0

#include "MSBVARcpp.h"

class SSobj
{
 private:

  // Private member data 
  int _T;           // number of observations
  int _h;           // number of regimes
  int _size;        // size of state-space (_T*_h)
  int _nint;        // int vals stored per regime [ceil(_T/16)] 
  int _st0;         // previous state 
  UINT **_int;  // state-space as 2D array _int[_h][_nint] 

public:
  
  int **_tmap;      // map of regime transitions

  SSobj(){_T=0;_h=0;_nint=0;_size=0;};    // Empty constructor
  SSobj(int, int);                        // SSobj ss(h,T);
  SSobj(SEXP);                            // Encode SS from R object  
  SSobj(UINT**,int**,int,int);        // Populate SS from arrays
 
  // Be a good steward of shared resources... pack out what you pack in 
  ~SSobj(){Free(_int); Free(_tmap);} 

  void encodeSS(SEXP);         // Encode SS as integer representation
  void SSinit(int, int);       // Initialize SSobj private member data 
  void SSbuild(UINT**, int**, int, int);
  void SSpopulate(UINT**, int**, int, int);
  void SSgenerate(Matrix &, Matrix &);
  void SSclean();
  void normalizeSS();

  // Accessor functions returning desired SS information   
  int SS_h();
  int SS_T();
  int SS_size();
  int SS_nint();
  int* SSinfo();

  void setregime(int, int); 
  int getregime(int); 
  SEXP getregimeR(int);

  int** getSSmap();
  SEXP getSSmapR();

  SEXP SSslim();          // R obj of encoded SS
  SEXP SSfull();          // R obj of full SS
};


class SSlist
{
 private:
  
  class SSobj* _SS;       // Array of _n SSobj objects 

 public:
  int _ndraws;            // Number of draws
  int _h;                 // Number of regimes
  int _T;                 // Number of observations
  
  // Constructors/Destructors

  SSlist(int); 
  SSlist(SEXP);
  //  ~SSlist(){Free(_SS);}

  // Accessor functions 

  void SSlistinit(int);
  void SSlistR(SEXP);
  void SSlistclean();
  void setSSobj(const class SSobj&, int);
  SSobj getSSobj(int);
  int** SSsum();
  double** SSmean(const int);
  double** SSvar(const int);
 
};

#endif

// Body file:  SS.cpp
