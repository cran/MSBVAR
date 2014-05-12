#include <R.h>
#include <Rmath.h>
#include "math.h"
#include "MSBVARcpp.h"

double getvlog(const Matrix &W, const Matrix &T, const ColumnVector &bf,
	       double cterm, int df, double tol)
{
  Matrix Vk, Vtr; ColumnVector gbeta; double lgbeta, ldetVtr, xbfVtri;
  Vk=T.i()*W; gbeta=Vk.i()*bf; QRD Vtrqr(Vk.t(),tol); Vtr=Vtrqr.R();
  lgbeta=log(fabs(gbeta(1))); //Rprintf("log(abs(gbeta[1])): %f\n", lgbeta);
  ldetVtr=log(fabs(Vtr.Determinant())); //Rprintf("ldetVtr: %f\n",ldetVtr);
  xbfVtri=(bf.t()*Vtr.i()*(bf.t()*Vtr.i()).t()).AsScalar(); //Rprintf("xbfVtri: %f\n",xbfVtri);
  return cterm-ldetVtr+df*lgbeta-0.5*df*xbfVtri;
}

//extern "C"
SEXP log_marg_A0k(SEXP WpostR, SEXP A0R, SEXP N2R, SEXP consttermR,
			     SEXP bfR, SEXP UTR, SEXP TinvR, SEXP dfR, SEXP n0R)
{
  int i, j, *dTi, db, m, N2, df, len;
  N2=INTEGER(N2R)[0]; df=INTEGER(dfR)[0]; m=INTEGER(coerceVector(listElt(WpostR,"m"),INTSXP))[0];

  double *pbfi, *lpa0, *cterm, lN2, tol, maxvlog, lqlog;
  lN2=log((double)N2); tol=1E-12; cterm=REAL(consttermR); //Rprintf("m: %d\nN2: %d\ndf: %d\n",m,N2,df);

  // Initialize Tinv/b.free/vlog variables
  SEXP Ti, bfi; Matrix Tinv; ColumnVector bfree, vlog(N2), qlog;

  // Initialize Wlist/W/Wmat objects and populate Wlist from WpostR
  Wlist Wall(WpostR,N2); Wobj W; Matrix Wmat;

  // Initialize SEXP/ptr to store/access log marginal A0k values
  SEXP lpa0yao; PROTECT(lpa0yao=allocVector(REALSXP,m)); lpa0=REAL(lpa0yao);
  for(i=0;i<m-1;i++){
    PROTECT(Ti=VECTOR_ELT(TinvR,i));
    dTi=getdims(Ti); Tinv=R2Cmat(Ti,dTi[0],dTi[1]);
    UNPROTECT(1); //Rprintf("Tinv[[%d]](%dx%d) initialized\n",i,dTi[0],dTi[1]);

    PROTECT(bfi=VECTOR_ELT(bfR,i));
    db=length(bfi); bfree.ReSize(db); pbfi=REAL(bfi); bfree<<pbfi;
    UNPROTECT(1); //Rprintf("bfree[[%d]](%d) initialized\n",i,db);

    for(j=1;j<=N2;j++){
      Wall.getWobj(W,j); Wmat=W.getWelt(i+1); W.clear();
      vlog(j)=getvlog(Wmat,Tinv,bfree,cterm[i],df,tol); //Rprintf("vlog(%d): %f\n",j,vlog(j));
    }

    // Modified harmonic mean of the max
    maxvlog=vlog.Maximum(); qlog=vlog-maxvlog; len=qlog.Storage();
    lqlog=0; for(j=1;j<=len;j++) lqlog+=exp(qlog(j)); lqlog=log(lqlog); // log(sum(exp(qlog)))
    lpa0[i]=maxvlog-lN2+lqlog; //Rprintf("lpa0[%d] = %f\n", i, lpa0[i]);
  }

  // Computations for last column
  PROTECT(Ti=VECTOR_ELT(TinvR,m-1));
  dTi=getdims(Ti); Tinv=R2Cmat(Ti,dTi[0],dTi[1]);
  UNPROTECT(1); //Rprintf("Tinv[[%d]](%dx%d) initialized\n",i,dTi[0],dTi[1]);

  PROTECT(bfi=VECTOR_ELT(bfR,m-1));
  pbfi=REAL(bfi); bfree.ReSize(length(bfi)); bfree<<pbfi;
  UNPROTECT(1); //Rprintf("bfree[[%d]](%d) initialized\n",i,db);

  UTobj UT(UTR); Matrix A0=R2Cmat(A0R,m,m);
  A0=drawA0cpp(A0,UT,df,INTEGER(n0R),W); Wmat=W.getWelt(m);
  lpa0[m-1]=getvlog(Wmat,Tinv,bfree,cterm[m-1],df,tol); //Rprintf("lpa0[%d] = %f\n",m-1,lpa0[m-1]);

  // Return R object lpa0yao
  UNPROTECT(1); return lpa0yao;
}
