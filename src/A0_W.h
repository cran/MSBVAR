#ifndef A0_W_H
#define A0_W_H

#include "MSBVARcpp.h"

/////////////////////////////////////////
// A0obj: Class interface for A0 draws //
/////////////////////////////////////////

class A0obj{
  // Private data members for A0 objects 
  int _m;                // # of eqn in sys
  int _nident;           // # of identified relationships
  int *_str;             // A0 structure
  int *_fpidx;           // Indices of free parameters 
  int *_xidx;            // Indices of restrictions 
  double *_val;          // A0 values
  
  // Private data members to track draws 
  int _ndraw;            // # of A0s stored
  int _cdraw;            // current A0 index

 public:
  A0obj(); ~A0obj();
  A0obj(SEXP); 

  // A0 constructor with structure and number of draws
  A0obj(const Matrix&, const int); 

  // A0 utility functions 
  void setA0(const Matrix&,const int); // Set A0 at supplied index
  void setA0c(const Matrix&);           // Set A0 at current index

  // A0 accessor functions
  int fpidx(const int) const;          // Free parameter idxs
  int xidx(const int) const;           // Restriction idxs
  int numfp() const;                   // # of free params
  ReturnMatrix structure() const;      // Get copy of A0 structure
  ReturnMatrix getA0(const int) const; // Get A0 at supplied index
  ReturnMatrix getA0c() const;          // Get A0 at current index
  double* getA0vals() const;           // Get copy of A0 value array
  SEXP toR() const;                    // Build R object from member data 
};

////////////////////////////////////
// Wobj: class to store W objects //
////////////////////////////////////

class Wobj{
  // Private data members for W objects 
  int _m;                // # of eqn in sys
  int *_idx;             // W matrix indicies 
  int _nval;
  double *_val;          // W vals flattened

 public:

  // W constructors/destructor
  Wobj(); 
  Wobj(SEXP);                               // Wobj(W) 
  Wobj(int, double*, int*);                 // Wobj(Wvals,m,Widx)
  ~Wobj();
  
  // W assignment functions
  void setW(SEXP);                          // Set W from R obj 
  void setWflat(int, double *, int *);
  void setWelt(const Matrix &, const int);        // Set W[[i]] 

  // W accessor functions
  SEXP getW() const;                        // Get W as R obj 
  double* getWvals() const;               // Get W as flattened array 
  int* getWindex() const;                   // Get W indicies 
  ReturnMatrix getWelt(int) const;    // Get W[[i]] as Matrix

  // Clear W member data 
  void clear();
};

///////////////////////////////////////////////////////////
// Wlist: class for storing/returning lists of W objects //
///////////////////////////////////////////////////////////

class Wlist : public Wobj {

  int _n;       // # draws
  int _m;       // # eq
  int _cdraw;   // current index
  int _nidx;    // # _idx elements
  int *_idx;    // _idx[_n*_m] array of indicies
  int _nval;    // # _val elements
  double *_val; // flat vector of W obj values

 public:
  Wlist(); ~Wlist();                     // Empty constructor
  Wlist(const int, const int);           // Wlist(num_draws, num_eq)
  Wlist(SEXP, const int);                // Wlist from R obj W.posterior
  SEXP toR() const;                      // Return formatted R obj

  int num_eq() const;                    // Return number of equations (_m)
  int num_draws() const;                 // Return number of draws (_n)

  void getWobj(Wobj&, int) const;         // Get W at index
  void setWobj(Wobj&, const int);        // Insert W into Wlist at index
};

/////////////////////////////////////////////////
// UTobj: class for gibbs.setup list object UT // 
/////////////////////////////////////////////////

class UTobj {
  int _m;                // num elements in UT 
  int _nval;             // num vals in UT
  int *_dim;        // dimensions
  double *_val;          // values in UT obj

 public:
  UTobj(); ~UTobj(); UTobj(SEXP); 

  void setUTelt(Matrix &, const int); 
  ReturnMatrix getUTelt(const int) const; 
  SEXP toR() const; 
}; 

#endif

////////////////////
// body: A0_W.cpp //
////////////////////
