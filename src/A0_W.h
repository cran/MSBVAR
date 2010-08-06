#ifndef A0_W_H
#define A0_W_H

#include "MSBVARcpp.h"

/////////////////////////////////////////
// A0obj: Class interface for A0 draws //
/////////////////////////////////////////

class A0obj{
  // Private data members for A0 objects 
  int eqnCount;//_m;                // # of eqn in sys
  int idRltnshpCount;//_nident;           // # of identified relationships
  int *A0Struct;//*_str;             // A0 structure
  int *freeParamIDX;//*_fpidx;           // Indices of free parameters 
  int *restrictionIDX;//*_xidx;            // Indices of restrictions 
  double *A0Values;//*_val;          // A0 values
  
  // Private data members to track draws 
  int numDraws;//_ndraw;            // # of A0s stored
  int currentDraw;//_cdraw;            // current A0 index

 public:
  A0obj(); //~A0obj();
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
  int eqnCount;                // # of eqn in sys
  int *WMatIDX;//*_idx;             // W matrix indicies 
  int numVals;//_nval;
  double *flatWVals;//*_val;          // W vals flattened

 public:

  // W constructors/destructor
  Wobj(); 
  Wobj(SEXP);                               // Wobj(W) 
  Wobj(int, double*, int*);                 // Wobj(Wvals,m,Widx)
  //~Wobj();
  
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

  int numDraws;//_n;       // # draws
  int eqnCount;       // # eq
  int currentDraw;//_cdraw;   // current index
  int IDXcount;//_nidx;    // # _idx elements
  int *arrayIDX;//*_idx;    // _idx[_n*_m] array of indicies
  int numVals;//_nval;    // # _val elements
  double *flatWValList;//*_val; // flat vector of W obj values

 public:
  Wlist(); //~Wlist();                     // Empty constructor
  Wlist(const int, const int);           // Wlist(num_draws, num_eq)
  Wlist(SEXP, int);                      // Wlist from R obj W.posterior
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
  int eleCount;                // num elements in UT 
  int numVals;//_nval;             // num vals in UT
  int *dimensionList;//_dim;        // dimensions
  double *objVals;//*_val;          // values in UT obj

 public:
  UTobj(); //~UTobj(); 
  UTobj(SEXP); 

  void setUTelt(Matrix &, const int); 
  ReturnMatrix getUTelt(int) const; 
  SEXP toR() const; 
}; 

#endif

////////////////////
// body: A0_W.cpp //
////////////////////
