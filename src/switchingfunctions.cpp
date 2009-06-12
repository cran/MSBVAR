#include "MSBVARcpp.h"

////////////////////////////////////////////////////////////////////
// 2008-11-29 : PTB updated bingen() and SSdraw so that it correctly does
// the filtering-smoothing-sampling steps.
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
// SSdraw() : Draws State-Spaces for MSBVAR models in C
////////////////////////////////////////////////////////////////////

extern "C" SEXP SSdraw(SEXP SSe, SEXP Xik, SEXP QR, SEXP SQR, SEXP TTR, SEXP hR, SEXP mR)
{
  // Allocations / declarations
  int TT=INTEGER(TTR)[0], h=INTEGER(hR)[0], m=INTEGER(mR)[0];
  Matrix Q = R2Cmat(QR, h, h);

  Matrix fp = BHLK(SSe, Xik, Q, SQR, TT, h, m);   // filtering step

  //  Rprintf("Data have been filtered\n");
  SSobj ss(fp.Ncols(), fp.Nrows());  
  SEXP ssout, out, names, tmap;; 
  int *sstmp;

  ss.SSgenerate(fp,Q);    // smoothing and sampling step
  //  Rprintf("State-space has been simulated!\n");

  // format the output of the 0-1 matrix
  PROTECT(ssout=allocMatrix(INTSXP, TT, h));
  sstmp = INTEGER(ssout);
  int i, j, tmp=0; 
  for(i=0; i<TT; i++) { 
    for(j=0; j<h; j++) {
      tmp=ss.getregime(i)-1;  // Get the regime
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
  //sp=sp/sp.Sum();

  int state=-1;
  //  Rprintf("Predicted Probs\n");
  // printMatrix(sp);
  for(i=1; i<=h-1; i++) 
    {
      Real pp = sp.Row(i).AsScalar()/sp.Rows(i,h).Sum(); // cumsum for state i
      double u = runif(0,1);
      //      Rprintf("Value of pp %f\n", pp);
      if(pp>u) state=i-1; break;    // sampled value
    }
  if(state>-1) {return state; }else{ return h-1;}
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
    x=u.Column(i).AsMatrix(TT,m);  xcov=sigma.Column(i).AsMatrix(m,m); 
    pYSt.Column(i) = dmvnorm(x,xcov); 
  }
  
  // Combat the threat of underflow
  double *pYStptr = pYSt.Store(), pYStmin=pow(10,-30); 
  for(i=0;i<pYSt.Storage();i++) if(pYStptr[i]<pYStmin) pYStptr[i]=pYStmin;


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
  pfilter.Release(); return pfilter.ForReturn();
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

// ////////////////////////////////////////////////////////////////////
// // Repeat each element of x 'each' times 
// // -------------------------------------
// // x: Vector of values to repeat
// // each: num of times to repeat each value in x
// ////////////////////////////////////////////////////////////////////

// ReturnMatrix rep(const ColumnVector& x, int each)
// {
//   int i, j, nrows=x.Storage(); ColumnVector repx(nrows*each); 
//   for(i=1;i<=nrows;i++) for(j=1;j<=each;j++) repx((i-1)*nrows+j)=x(i);
//   repx.Release(); return repx.ForReturn(); 
// }

////////////////////////////////////////////////////////////////////
// Multivariate normal density function
////////////////////////////////////////////////////////////////////

ReturnMatrix dmvnorm(Matrix &x, Matrix &sigma)
{
  int i, n=x.Nrows();  ColumnVector dist(n); 
  Matrix xsigmai=x*sigma.i();  xsigmai=SP(xsigmai,x); 
  double dtmp = x.Ncols()*LOG2PI+sigma.LogDeterminant().LogValue();
  for(i=1;i<=n;i++) dist(i)=exp(-(dtmp+xsigmai.Row(i).Sum())/2.0); 
  dist.Release(); return dist.ForReturn();
}


// ////////////////////////////////////////////////////////////////////
// // SSsumR
// ////////////////////////////////////////////////////////////////////
// SEXP SSsumR(SEXP SSsampleR)
// {
//   int i, j;  SSlist SSsample(SSsampleR); SEXP out; 
//   PROTECT(out=allocMatrix(INTSXP, SSsample._T, SSsample._h));
//    int *pout=INTEGER(out), **tmp=SSsample.SSsum(); 
//    for(i=0;i<SSsample._h;i++) {
//      for(j=0;j<SSsample._T;j++) {
//        pout[i*SSsample._T+j]=tmp[i][j];
//      }}

//   UNPROTECT(1);
//   Free(tmp);
//   return out;
// }

// ////////////////////////////////////////////////////////////////////
// // SSmeanR
// ////////////////////////////////////////////////////////////////////

// SEXP SSmeanR(SEXP SSsampleR, SEXP method)
// {
//   int i, j, m;  double **means, *pout; 
//   SSlist SSsample(SSsampleR); SEXP out;  
  
//   m=INTEGER(method)[0]; 
//   Rprintf("m=%d\t_T=%d\t_h=%d\n",m,SSsample._T,SSsample._h);
//   means=SSsample.SSmean(m); 

//   switch(m){
//   default:
//   case 1: {
//     PROTECT(out=allocMatrix(REALSXP, SSsample._T, SSsample._h));
//     pout=REAL(out); 
//     for(i=0;i<SSsample._h;i++) 
//       for(j=0;j<SSsample._T;j++) 
// 	pout[i*SSsample._T+j]=means[i][j];
//     break;
//   }

// //   case 2: {
// //     PROTECT(out=C2Rdoublemat(means));
// //     break;
// //   }
//   }

//   Free(means);  UNPROTECT(1);  
//   return out;
// }

// SEXP SSvarR(SEXP SSsampleR, SEXP method)
// {
//   int i, j, m, h; 
//   SSlist SSsample(SSsampleR); 
//   m=INTEGER(method)[0]; h=SSsample._h;
//   SEXP out;  PROTECT(out=allocVector(REALSXP, h)); 
//   double **var=SSsample.SSvar(m), *pout=REAL(out); 
//   for(i=0;i<SSsample._h;i++) for(j=0;j<SSsample._T;j++) pout[i*SSsample._T+j]=var[i][j];
//   UNPROTECT(1); return out;
// }
