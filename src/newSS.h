#ifndef newSS_H
#define newSS_H
//#define UINT uint32_t
//#define NBITS 30
//#define NBITSF 30.0

#include "MSBVARcpp.h"

class newSSobj
{
  private:

  // Private member data 
  int T;           	// number of observations
  int H;           	// number of regimes
  int Size;        	// size of state-space (observations*regimes)
  int prevState;   	// previous state 
  int **StateSpace;  	// state-space as 2D array int[regimes][observations] 

  public:
  
  int **TransMap;      // map of regime transitions

  newSSobj(){T=0;H=0;Size=0;};		// Empty constructor
  newSSobj(int, int);          		// SSobj ss(reg,obs);
  //newSSobj(int**,int**,int,int);        // Populate SS from arrays
 
  //Deconstruction handled by R
  
  // Accessor functions returning desired SS information
  int countSSReg(); 	//Get H
  int countSSObs(); 	//Get T
  int getSize();	//Get H*T
  bool hasEachObs();	//Returns whether or not an observation has been made in each state.
  SEXP getSSMapR(); 	//R version of SS map
  SEXP getTransMapR();	//R version of Trans map
  //SEXP getSEXP(SEXP);	//Formats SS/Trans for R return as a single array.

  //Constructor functions for creating an SS object
  void Build(int, int); //Construction method given regime, observations count
  //void Build(int, int, int**, int**); 	//Fully constructive method

  // Modification functions. setRegime is for transitions on the state space
  // getRegime returns a regime at an int time
  void setFirst(int, int);
  void setRegime(int, int); 
  void reset();		//Zeroes out SS and Trans.
  int getRegime(int); 	//Returns a regime at a given time.
};
#endif
// Body file:  newSS.cpp
