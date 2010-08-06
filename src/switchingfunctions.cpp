#include "MSBVARcpp.h"
#include <Rdefines.h>

////////////////////////////////////////////////////////////////////
// 2008-11-29 : PTB updated bingen() and SSdraw so that it correctly does
// the filtering-smoothing-sampling steps.
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
// SSdraw() : Draws State-Spaces for MSBVAR models in C
////////////////////////////////////////////////////////////////////

//New generation function to perform FFBS 
//to replace SSGenerate
extern "C" SEXP newSSdraw(SEXP SSe, SEXP Xik, SEXP QR, SEXP SQR, SEXP TTR, SEXP hR, SEXP mR)
{
  GetRNGstate();	//Necessary to avoid inconsistency

  // Allocations / declarations
  int TT=INTEGER(TTR)[0], h=INTEGER(hR)[0], m=INTEGER(mR)[0];
  Matrix Q = R2Cmat(QR, h, h);

  Matrix fp = BHLK(SSe, Xik, Q, SQR, TT, h, m);   // filtering step

  //Allocate SS object and associated variables
  //out is the return object, names are $SS and such
  //ssout holds SS, tmap holds the trans map.
  SEXP ssout, out, names, tmap;
  //T is observations, i is an index variable
  //h2 is a regime storage variable, state holds current iterative state
  int T, i, h2, state;
  Matrix tmp;	//Used for 'coinflip'

  newSSobj ss(fp.Ncols(), fp.Nrows());	//Declare our ss object, manipulated in do-while
  do
  {
	ss.reset();		//Zero out ss and trans
	T=ss.countSSObs();	//T holds observation count
	tmp = fp.Row(T);	//Allocate tmp matrix
	h2=Q.Nrows();		//h2 holds regime count
	state=-1;		//state defaults to -1 for later if statement
  	
	//Pre-bingen run
  	for(i=1; i<=h2-1; i++) //Iterate across regimes
  	{
	    //"Coinflip" happens here
  	    Real pp = tmp.Column(i).AsScalar()/tmp.Columns(i,h2).Sum();
  	    double u = runif(0,1);
  	    if(pp<=u) 
  	    {
		state=i-1; 
	      	break; 
	    }   
	}
	if(state>-1)	//If we found a state
	  ss.setFirst(T, state); //setFirst used to avoid extra transition
	else		//Otherwise, default to final regime
	  ss.setFirst(T, h2-1); 
	  
	// Now back-sample from T-1:1
	for(int t=T-1; t>0; t--)
	{ 
	   // Get the filter probability and smooth it and sample it
	   tmp=fp.Row(t); 
	   ss.setRegime(t-1, bingen(tmp, Q, ss.getRegime(t)-1)); 
	}
  }
  while(!ss.hasEachObs());	//Do it again until all states have an observation.

  PROTECT(tmap=ss.getTransMapR());  	// get the transition matrix
  PROTECT(ssout=ss.getSSMapR());	// and the state space matrix
  PROTECT(out=allocVector(VECSXP, 2));	// Prepare a vector to receive the matrices
  PROTECT(names=allocVector(STRSXP, 2));// And another one to hold their names

  SET_VECTOR_ELT(out, 0, ssout);	//Give out the state space matrix
  SET_STRING_ELT(names, 0, mkChar("SS"));//And name it $SS
  SET_VECTOR_ELT(out, 1, tmap);		//Then give it the trans map
  SET_STRING_ELT(names, 1, mkChar("transitions"));//And name it $transitions

  setAttrib(out, R_NamesSymbol, names);//Tie the names and out together

  PutRNGstate();	//Put the RNG state(DO THIS BEFORE YOU UNPROTECT)

  UNPROTECT(4);

  //In case the garbage collector gives us more trouble,
  //The following code can be used with newSS's getSEXP
  //Function to pile SS and Trans into a large integer array.
  /*int reg = ss.countSSReg();
  int obs = ss.countSSObs();
  PROTECT(out = NEW_INTEGER(reg*obs*reg*reg));
  out = ss.getSEXP(out);
  PutRNGstate();
  UNPROTECT(1);*/

  return(out);
}


////////////////////////////////////////////////////////////////////
// Function to generate the regime of the state    
// space for one observation; returned as an integer.
// --------------------------------------------------
// p = h filtered probabilities of each regime for an observation         
// Q = h x h transition matrix
// st0 = previous regime index
////////////////////////////////////////////////////////////////////

int bingen(Matrix &p, Matrix &Q, int st0)
{
  int i=1, h=Q.Nrows(); 
  Matrix sp=SP(Q.Column(st0+1),p.t());  // one-step prediction t|t

  int state=-1;
  for(i=1; i<=h-1; i++) 
  {
      Real pp = sp.Row(i).AsScalar()/sp.Rows(i,h).Sum(); // cumsum for state i

      double u = runif(0,1);

      if(pp>u) 
      {
	state=i-1; 
      	break; 
      }   // sampled value
  }
  if(state>-1) 
    return state;
  else
    return h-1;
}

////////////////////////////////////////////////////////////////////
// BHLK Filter 
// -----------
// uR[TT*m x h] = residuals for reg model in each regime 
// sigmaR[m*m x h] = cov of residuals
// Q[h x h] = transition probability matrix
// SQR[h x 1] = steady state values of the chain
// TT = number of obs
// h = number of states
// m = number of variables
////////////////////////////////////////////////////////////////////

ReturnMatrix BHLK(SEXP uR, SEXP sigmaR, Matrix& Q, SEXP SQR, int TT, int h, int m)
{
  // Setup Constants
  int  i, t;

  // Setup input args and storage
  Matrix u=R2Cmat(uR,TT*m,h), sigma=R2Cmat(sigmaR,m*m,h), SQ=R2Cmat(SQR, h, 1);
  Matrix pYSt(TT,h), pfilter(h,TT), nptop(h,h), utmp(TT,m), sigtmp(m,m);
  //ColumnVector steady(h); 
  Matrix x(TT,m), xcov(m,m);

  // Compute the likelihood for each obs. in each regime 
  for(i=1; i<=h; i++)
  {
    x=u.Column(i).AsMatrix(TT,m);  
    xcov=sigma.Column(i).AsMatrix(m,m); 
    pYSt.Column(i) = dmvnorm(x,xcov); 
  }
  
  // Combat the threat of underflow
  double *pYStptr = pYSt.Store(), pYStmin=pow(10.0,-30); 
  for(i=0;i<pYSt.Storage();i++) 
    if(pYStptr[i]<pYStmin) 
      pYStptr[i]=pYStmin;


  // Set up past state
  ColumnVector steady = Q*SQ;

  // Here's the loop to compute the filtered probabilities -- note
  // that this is vectorized, unlike earlier versions
  for(t=1; t<=TT; t++)
  {
    ColumnVector num = SP(pYSt.Row(t).t(), steady);
    Matrix den = pYSt.Row(t) * steady;
    pfilter.Column(t) = num/den.AsScalar();
    ColumnVector steady = Q*(num/den.AsScalar());
  }

  pfilter = pfilter.t();
  return pfilter.ForReturn();
}

////////////////////////////////////////////////////////////////////
// BHLK filter function wrapper 
// ----------------------------
// 1) callable from R via BHLK.filter(u, Sigma, Q)
// 2) returns classed R object  
////////////////////////////////////////////////////////////////////

SEXP BHLKR(SEXP uR, SEXP SigmaR, SEXP QR, SEXP SQR, SEXP TR, SEXP hR, SEXP mR)
{ 
  int TT=INTEGER(TR)[0], h=INTEGER(hR)[0], m=INTEGER(mR)[0];
  //  Matrix u=R2Cmat(uR,TT*m,h);
  Matrix Q=R2Cmat(QR,h,h), filtered=BHLK(uR, SigmaR, Q, SQR, TT, h, m);
  return C2Rmat(filtered); 
}

////////////////////////////////////////////////////////////////////
// Multivariate normal density function
////////////////////////////////////////////////////////////////////

ReturnMatrix dmvnorm(Matrix &x, Matrix &sigma)
{
  int i, n=x.Nrows();  
  ColumnVector dist(n); 
  Matrix xsigmai=x*sigma.i();  
  xsigmai=SP(xsigmai,x); 
  double dtmp = x.Ncols()*LOG2PI+sigma.LogDeterminant().LogValue();
  for(i=1;i<=n;i++) 
    dist(i)=exp(-(dtmp+xsigmai.Row(i).Sum())/2.0); 
  return dist.ForReturn();
}


//Should be removed in the immediate future.
/*extern "C" SEXP SSdraw(SEXP SSe, SEXP Xik, SEXP QR, SEXP SQR, SEXP TTR, SEXP hR, SEXP mR)
{
  // Allocations / declarations
  int TT=INTEGER(TTR)[0], h=INTEGER(hR)[0], m=INTEGER(mR)[0];
  Matrix Q = R2Cmat(QR, h, h);

  Matrix fp = BHLK(SSe, Xik, Q, SQR, TT, h, m);   // filtering step

  //  Rprintf("Data have been filtered\n");
  SSobj ss(fp.Ncols(), fp.Nrows());  
  SEXP ssout, out, names, tmap;
  int *sstmp;

  ss.SSgenerate(fp,Q);    // smoothing and sampling step for the SS class
  
  // format the output of the 0-1 matrix
  PROTECT(ssout=allocMatrix(INTSXP, TT, h));
  sstmp = INTEGER(ssout);

  int i, j, tmp=0; 

  for(i=0; i<TT; i++) 
  { 
    for(j=0; j<h; j++) 
    {
      tmp=ss.getregime(i)-1;  // Get the regime
      //Rprintf("Regime :%d\n", tmp);
      
      sstmp[i + TT*j] = 0;    // Zero things out
    }
    sstmp[i + TT*tmp] = 1;    // Fill in the 1
  }

  PROTECT(tmap=ss.getSSmapR());  // get the transition matrix

  PROTECT(out=allocVector(VECSXP, 2));
  PROTECT(names=allocVector(STRSXP, 2));
  
  SET_VECTOR_ELT(out, 0, ssout);
  SET_STRING_ELT(names, 0, mkChar("SS"));
  SET_VECTOR_ELT(out, 1, tmap);
  SET_STRING_ELT(names, 1, mkChar("transitions"));

  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(4);

  return(out);
}*/
