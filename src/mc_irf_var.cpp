#include "MSBVARcpp.h"

//extern "C"
SEXP mc_irf_var(SEXP varobj, SEXP nsteps, SEXP draws)
{
  int m, p, dr=INTEGER(draws)[0], ns=INTEGER(nsteps)[0], T, df, i;
  SEXP AR, Y, Bhat, XR, prior, hstar, meanS, output;

  // Get # vars/lags/steps/draws/T/df
  PROTECT(AR = listElt(varobj, "ar.coefs"));
  PROTECT(Y = listElt(varobj, "Y"));
  m = INTEGER(getAttrib(AR, R_DimSymbol))[0]; //#vars
  p = INTEGER(getAttrib(AR, R_DimSymbol))[2]; //#lags
  T = nrows(Y); df = T - m*p - m - 1;
  UNPROTECT(2);

  // Put coefficients from varobj$Bhat in Bcoefs vector (m^2*p, 1)
  PROTECT(Bhat = coerceVector(listElt(varobj, "Bhat"), REALSXP));
  Matrix bcoefs = R2Cmat(Bhat, m*p, m);
  bcoefs = bcoefs.AsColumn();
  UNPROTECT(1);

  // Define X(T x m*p) subset of varobj$X and XXinv as solve(X'X)
  PROTECT(XR = coerceVector(listElt(varobj,"X"),REALSXP));
  Matrix X = R2Cmat(XR, T, m*p), XXinv;
  UNPROTECT(1);

  // Get the correct moment matrix
  PROTECT(prior = listElt(varobj,"prior"));
  if(!isNull(prior)){
    PROTECT(hstar = coerceVector(listElt(varobj,"hstar"),REALSXP));
    XXinv = R2Cmat(hstar, m*p, m*p).i();
    UNPROTECT(1);
  }
  else { XXinv = (X.t()*X).i(); }
  UNPROTECT(1);

  // Get the transpose of the Cholesky decomp of XXinv
  SymmetricMatrix XXinvSym; XXinvSym << XXinv;
  XXinv = Cholesky(XXinvSym);

  // Cholesky of covariance
  PROTECT(meanS = coerceVector(listElt(varobj,"mean.S"),REALSXP));
  SymmetricMatrix meanSSym; meanSSym << R2Cmat(meanS, m, m);
  Matrix Sigmat = Cholesky(meanSSym);
  UNPROTECT(1);

  // Matricies needed for the loop
  ColumnVector bvec; bvec=0.0;
  Matrix sqrtwish, impulse(dr,m*m*ns); impulse = 0.0;
  SymmetricMatrix sigmadraw; sigmadraw = 0.0;
  IdentityMatrix I(m);

  GetRNGstate();
  // Main Loop
  for (i=1; i<=dr; i++){
    // Wishart/Beta draws
    sigmadraw << Sigmat*(T*rwish(I,df).i())*Sigmat.t();
    sqrtwish = Cholesky(sigmadraw);
    bvec = bcoefs+KP(sqrtwish, XXinv)*rnorms(m*m*p);

    // IRF computation
    impulse.Row(i) = irf_var_from_beta(sqrtwish, bvec, ns).t();
    if (!(i%1000)){ Rprintf("Monte Carlo IRF Iteration = %d\n",i); }
  } // end main loop
  PutRNGstate();

  int dims[]={dr,ns,m*m};
  PROTECT(output = C2R3D(impulse,dims));
  setclass(output,"mc.irf.VAR");
  UNPROTECT(1);
  return output;
}

