#include "MSBVARcpp.h"

//extern "C"
SEXP mc_irf_bsvar(SEXP A0postR, SEXP nstepsR, SEXP N2R, SEXP mR,
			     SEXP pR, SEXP ncoefR, SEXP n0R, SEXP n0cumR,
			     SEXP XXinvR, SEXP Ui, SEXP PpostR, SEXP signlistR)
{
  int i, j, nsteps, N2, m, p, ncoef, pctct=0, *n0, *n0cum;
  nsteps=INTEGER(nstepsR)[0]; N2=INTEGER(N2R)[0];
  m=INTEGER(mR)[0]; p=INTEGER(pR)[0]; ncoef=INTEGER(ncoefR)[0];
  n0=INTEGER(n0R); n0cum=INTEGER(n0cumR);
  //  Rprintf("m = %d\np = %d\nncoef = %d\nnsteps = %d\nN2 = %d\n", m, p, ncoef, nsteps, N2);

  A0obj A0post(A0postR); Matrix A0(m,m), A0i(m,m); //   Rprintf("A0 variables assigned\n");
  Matrix XXinv=R2Cmat2(XXinvR), IRF(N2,m*m*nsteps); //   Rprintf("XXinv variable assigned\n");
  ColumnVector b; RowVector signvec=R2CRV(signlistR);
  Matrix signmat(m,m); for(i=1;i<=m;i++) signmat.Row(i)=signvec;
//   Rprintf("signmat(%dx%d):",signmat.Nrows(),signmat.Ncols()); printMatrix(signmat);

  GetRNGstate();
  for(i=0;i<N2;i++){
    // Get A0 object from A0.posterior
    A0=A0post.getA0(i); //     Rprintf("A0(%dx%d)\n",A0.Nrows(),A0.Ncols()); printMatrix(A0);

    // Schur product of A0 and the matrix of signs for renormalization
    A0=SP(A0,signmat);

    A0i=A0.i(); //     Rprintf("A0i(%dx%d)\n",A0i.Nrows(),A0i.Ncols()); printMatrix(A0i);
    b=a2b(A0,Ui,n0,n0cum); //     Rprintf("b(%d)\n\n",b.Storage());

    Matrix F(ncoef,m); F=0.0;
    for(j=0;j<m;j++)
      F.Column(j+1)=R2Cmat2(VECTOR_ELT(PpostR,j))*b.Rows(n0cum[j]+1,n0cum[j+1]);

    F+=XXinv*rnorms_mat(ncoef,m);
    F=(F.Rows(1,m*p)*A0i).t(); //printMatrix(F);
    b=F.AsRow().t();      //printCVector(b);

    IRF.Row(i+1)=irf_var_from_beta(A0i.t(), b, nsteps).t(); //     Rprintf("IRF.Row(i+1)=irf_var_from_beta(A0i.t(), b(%d), nsteps); \n",b.Storage());
    if(i%(N2/10)==0) Rprintf("Monte Carlo IRF %d %% Complete\n",++pctct*10);
  }
  PutRNGstate();

  SEXP out; int dIRF[]={N2,nsteps,m*m};
  PROTECT(out=C2R3D(IRF,dIRF));
  setclass(out,"mc.irf.BSVAR");
  UNPROTECT(1); return out;
}
