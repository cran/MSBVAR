#include "newSS.h"
#include "MSBVARcpp.h"
#include "Rdefines.h"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
////                                                  ////
////  Functions for class SSobj: an ANSI C compatible ////
////  implementation of state-spaces object storage   ////
////  and processing.                                 ////
////                                                  ////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

// Constructors that populate integer state space matrix.
//Effectively required wrappers for the Build methods.
/*newSSobj::newSSobj(int **ss, int **tmap, int reg, int obs)
{
  Build(reg, obs, ss, tmap);
}*/

newSSobj::newSSobj(int reg, int obs)
{ 
  Build(reg, obs);
}

//Construction assistance methods
void newSSobj::Build(int reg, int obs)
{
  H = reg; //Regimes
  T = obs; //Observations
  Size = reg*obs; 
  prevState = 0; //prevState defaults to zero

  //Allocation happens here, never use anything but R_alloc
  //Use of other allocation functions yields segfaults
  StateSpace = (int **)R_alloc(H, sizeof(int*));
  TransMap = (int **)R_alloc(H, sizeof(int*)); 

  // Initialize arrays to zero
  // I'm not sure this is necessary, R_alloc should initialize to all bits to 0.
  // Could be worth removing to see if anything breaks in the future.
  // I mean, it's a lot of potential processing reduction.
  int i, j;
  for(i=0;i<H;i++)
  { 
    StateSpace[i] = (int *)R_alloc(T, sizeof(int)); 
    for(j=0;j<T;j++) //State Space zeroing
      StateSpace[i][j] = 0;
    TransMap[i] = (int *)R_alloc(H, sizeof(int));  
    for(j=0;j<H;j++) //Transition map zeroing
      TransMap[i][j] = 0;
  } 
}

//This particular constructive method is currently unused.
//It's generally just as easy to build the state and trans
//Matrices inside an int,int constructed newSS object.
//Uncomment if we need to be able to do this for some reason.
/*void newSSobj::Build(int reg, int obs, int** state, int** trans)
{
  H = reg; //Regimes
  T = obs; //Observations
  Size = reg*obs;
  prevState = 0;
  int i, j;
  //Again, don't use anything but R_alloc, dragons etc
  StateSpace = (int **)R_alloc(H, sizeof(int*));
  TransMap = (int **)R_alloc(H, sizeof(int*));
  
  //Same for-loop as the 0-out above, just uses state and trans instead of 0.
  for(i=0;i<H;i++)
  {
    StateSpace[i] = (int *)R_alloc(T, sizeof(int));
    for(j=0;j<T;j++) 
	StateSpace[i][j]=state[i][j];
    TransMap[i] = (int *)R_alloc(H, sizeof(int));  
    for(j=0;j<H;j++) 
	TransMap[i][j]=trans[j][i];
  }
}*/

//Return regime count
int newSSobj::countSSReg()
{
  return H;
}

//Return observation count
int newSSobj::countSSObs()
{
  return T;
}

//Return state space size(h*T)
int newSSobj::getSize()
{
  return Size;
}

//Could feasibly consolidate the two below into one function
//That accepts a couple variables. Not necessary, though.
//R representation of the Trans map.
SEXP newSSobj::getTransMapR()
{
  SEXP out; 
  PROTECT(out=allocMatrix(INTSXP,H,H));
  int i, j, *pout=INTEGER(out);
  for(i=0;i<H;i++) 
	for(j=0;j<H;j++) 
		pout[j+i*H]=TransMap[i][j]; 
  UNPROTECT(1); 
  return out;
}

//R representation of the SS map
SEXP newSSobj::getSSMapR()
{
  SEXP out; 
  PROTECT(out=allocMatrix(INTSXP,T,H));
  int i, j, *pout=INTEGER(out);
  for(i=0;i<H;i++) 
	for(j=0;j<T;j++) 
		pout[j+i*T]=StateSpace[i][j]; 
  UNPROTECT(1); 
  return out;
}

//Special setRegime for initialization prior to bingen.
//Using normal setRegime instead of setFirst causes a
//"wandering transition", an extra transition. That's bad.
//Could be accomplished just as easily by adding a bool to
//setRegime, but this works too. If things get messy, do that.
void newSSobj::setFirst(int t, int st)
{
  StateSpace[st][t] = 1;	//Set regime at time t
  prevState = st;		//Update previous regime
}

// Set the bit for the regime
void newSSobj::setRegime(int t, int st)
{
  StateSpace[st][t] = 1;   	// Set regime at time t 
  TransMap[prevState][st]++;    // Update transition map 
  prevState=st;                 // Update previous regime 
}

//Zeroes out the state space and transition matrices.
void newSSobj::reset()
{
  int i, j;
  for(i=0;i<H;i++)
  { 
    for(j=0;j<T;j++) //State Space zeroing
      StateSpace[i][j] = 0;
    for(j=0;j<H;j++) //Transition map zeroing
      TransMap[i][j] = 0;
  } 
}

//Determines whether an observation exists for each state
//In the state space matrix.
bool newSSobj::hasEachObs()
{
	int tmp=0;	//Sum storage
	for(int i=0; i<H; i++)	//HxH loop through the trans map
	{
		for(int j=0;j<H;j++)	//For each row in the map
		{
			tmp += TransMap[i][j];	//Add the numbers together
		}
		if(tmp<1)	//If there is not at least one observation,
			return false;
		tmp=0;	//Reset your sum, loop.
	}//If you don't return false,
	return true; 
}

// Return the regime at time T as an integer
int newSSobj::getRegime(int t)
{
  int i;
  for(i=0; i<H-1; i++)
  {
    if(StateSpace[i][t] == 1)	//Does this regime transfer at t?
      return i+1;		//Yes, return the index immediately after it.
  }
  return H; //Otherwise, return the last possible regime's index.
}
//R representation of SS and Trans maps as a single
//Integer array. Useful for avoiding garbage collection
//Shenanigans if it becomes a problem(it shouldn't).
//'out' should be allocated by NEW_INTEGER for this.
//If it's not, bugs.
/*SEXP newSSobj::getSEXP(SEXP out)
{
	SEXP ret = out;
	int i,j;
	int index = 0;
	for(i=0;i<T;i++);
		for(j=0;j<H;j++)
		{
			INTEGER_POINTER(ret)[index] = StateSpace[i][j];
			index++;
		}
	for(i=0;i<T;i++)
		for(j=0;j<T;j++)
		{
			INTEGER_POINTER(ret)[index] = TransMap[i][j];
		}
	return ret;	
}*/
