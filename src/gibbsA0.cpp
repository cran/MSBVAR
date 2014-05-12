#include "MSBVARcpp.h"
#include "A0_W.h"

//extern "C"
SEXP gibbsA0(SEXP varobj, SEXP N1R, SEXP N2R, SEXP thinR,
			SEXP method, SEXP UTR)
{
  int i, j, k, m, df, norm_method=INTEGER(method)[0], pctct=0;
  int N1=INTEGER(N1R)[0], N2=INTEGER(N2R)[0], thin=INTEGER(thinR)[0], switchct=0;
  SEXP n0R, dfR, mode, identR;
  double *pA0;
  Matrix A0ml, A0gbs, ident;

//   Rprintf("N1:\t%d\nN2:\t%d\nthin:\t%d\nmethod:\t%d\n", N1, N2, thin, norm_method);

  // Initialize varobj$A0.mode/m/df/n0/ident
  PROTECT(mode=listElt(varobj,"A0.mode"));
  m=getdims(mode)[0];
  A0ml=R2Cmat(mode,m,m);
  UNPROTECT(1);

  PROTECT(dfR=coerceVector(listElt(varobj,"df"),INTSXP));
  df=INTEGER(dfR)[0];
  UNPROTECT(1);

  PROTECT(identR=listElt(varobj,"ident"));
  ident=R2Cmat(identR,m,m);
  double *pid=ident.Store();
  UNPROTECT(1);

  PROTECT(n0R=listElt(varobj,"n0"));
  int len=length(n0R), *pn0=INTEGER(n0R), *n0=(int*)R_alloc(len,sizeof(int));
  for(i=0;i<len;i++)
	n0[i]=pn0[i];
  UNPROTECT(1);

  // Normalize A0's mode/Set A0gbs
  A0ml=norm_svar(A0ml,A0ml,norm_method,&switchct);
  A0gbs=A0ml;

  // Initialize storage/loop variables
  A0obj A0all(ident, N2);
  Wlist Wall(N2,m);
  Wobj W;
  UTobj UT(UTR);

  // Find Matrix pointer locations of restrictions
  int numx=m*m-A0all.numfp(), *restx=(int*)R_alloc(numx,sizeof(int)), xct=0;
  for(i=0;i<m*m;i++)
	if(!pid[i])
		restx[xct++]=i;

  //////////////////
  // Burn-in Loop //
  //////////////////

  for(i=1;i<=N1;i++){
    // Draw posterior sample
    A0gbs=drawA0cpp(A0gbs, UT, df, n0, W);
    pA0=A0gbs.Store();

    // Clean up A0 to minimize roundoff error and normalize
    for(j=0;j<numx;j++)
	pA0[restx[j]]=0;
    A0gbs=norm_svar(A0gbs, A0ml, norm_method, &switchct);

    if(!(i%(N1/10)))
	Rprintf("Gibbs Burn-in %d %% Complete\n",++pctct*10);
  }

  // Reset counters
  switchct=0; pctct=0;

  ///////////////////
  // Sampling Loop //
  ///////////////////

  for(i=1;i<=N2;i++){
    // Take 'thin' number of draws
    for(j=1; j<=thin; j++)
    {
      // Draw posterior sample: Assign A0gbs explicitly, Wobj by reference
      A0gbs=drawA0cpp(A0gbs, UT, df, n0, W);
      pA0=A0gbs.Store();

      // Clean up A0 to minimize roundoff error and normalize
      for(k=0;k<numx;k++)
	pA0[restx[k]]=0;
      A0gbs=norm_svar(A0gbs, A0ml, norm_method, &switchct);
    }

    // Store A0/W in compact form
    A0all.setA0(A0gbs,i);
    Wall.setWobj(W,i);

    // Report iteration status and current A0's log determinant
    if(!(i%(N2/10))){
      Rprintf("Gibbs Sampling %d %% Complete (%d draws)\n", ++pctct*10, i);
      LogAndSign ld = A0gbs.LogDeterminant();
      Rprintf("A0 log-det \t = %f \n", ld.Sign()*ld.LogValue());
    }
  }

  // Create output list w/ A0.posterior, W.posterior, ident, thin, N2
  SEXP out, names;

  PROTECT(out=allocVector(VECSXP,5));
  PROTECT(names=allocVector(STRSXP,5));

  SET_VECTOR_ELT(out, 0, A0all.toR());
  SET_STRING_ELT(names, 0, mkChar("A0.posterior"));
  SET_VECTOR_ELT(out, 1, Wall.toR());
  SET_STRING_ELT(names, 1, mkChar("W.posterior"));
  SET_VECTOR_ELT(out, 2, C2Rmat(ident));
  SET_STRING_ELT(names, 2, mkChar("ident"));
  SET_VECTOR_ELT(out, 3, thinR);
  SET_STRING_ELT(names, 3, mkChar("thin"));
  SET_VECTOR_ELT(out, 4, N2R);
  SET_STRING_ELT(names, 4, mkChar("N2"));

  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(2);

  // Return output list
  return out;
}
