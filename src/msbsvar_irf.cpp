#include "MSBVARcpp.h"

SEXP msbsvar_irf(SEXP gibbs, SEXP msbsvar, SEXP nsteps)
{
  int i, k, n, N2, h, m, p, n0max, ns=INTEGER(nsteps)[0];
  int *db, *dF, *dxi, *dQ, N210pct, pctct=0;
  SEXP bR, FR, xiR, QR, Ui, IRFlist, IRFtmp;

//   Rprintf("ns = %d\n",ns); 
  
  // Get b, F, xi, Q, SS, dims from gibbs object
  PROTECT(bR = VECTOR_ELT(gibbs,0)); db=getdims(bR);
//   Rprintf("b(%d,%d)\n",db[0],db[1]); 
  PROTECT(FR = VECTOR_ELT(gibbs,1)); dF=getdims(FR);
//   Rprintf("F(%d,%d)\n",dF[0],dF[1]); 
  PROTECT(xiR= VECTOR_ELT(gibbs,2)); dxi=getdims(xiR);
//   Rprintf("xi(%d,%d)\n",dxi[0],dxi[1]); 
  PROTECT(QR = VECTOR_ELT(gibbs,3)); dQ=getdims(QR); UNPROTECT(1); 
//   Rprintf("Q(%d,%d)\n",dQ[0],dQ[1]); 
  
//   Rprintf("Gibbs Objects and Dimensions Assigned\n"); 
  
  // Reconstruct constants
  N2=db[0]; h=(int)sqrt((double)dQ[1]); n0max=db[1]/h; m=dxi[1]/h; p=((dF[1]/(h*m))-1)/m;
  N210pct=N2/10; 

//   Rprintf("N2=%d\nh=%d\nm=%d\np=%d\nn0max=%d\n",N2,h,m,p,n0max); 

  // Get Ui from msbsvar
  PROTECT(Ui=VECTOR_ELT(msbsvar,7));
 
  Matrix bsample=R2Cmat(bR,N2,n0max*h);
  Matrix Fsample=R2Cmat(FR,N2,m*(m*p+1)*h); 
  Matrix xisample=R2Cmat(xiR,N2,m*h);

  ColumnVector bk(n0max), Fk(m*(m*p+1)), bvec(m*m*p); bk=0.0; Fk=0.0; bvec=0.0; 
  DiagonalMatrix xik(m), sqrtxik(m); xik=0.0; sqrtxik=0.0; 
  Matrix Q(h,h), A0(m,m), A0i(m,m), fmat(m,m*p+1), sqrtwish, impulse(N2,m*m*ns); 
  double *pFk; int IRFdims[]={N2,ns,m*m};   

  PROTECT(IRFlist=allocVector(VECSXP,h)); 
  // Loop over regimes 
  for(k=1;k<=h;k++){
    
//     Rprintf("\n==========\nRegime %d\n==========\n",k);
    pctct=0;
    // Compute impulse responses for every draw of regime k
    for(n=1;n<=N2;n++){
//        Rprintf("\nDraw %d:\n",n); 

      // Get values for draw 'n', regime 'k' 
      bk=bsample.SubMatrix(n,n,(k-1)*n0max+1,k*n0max).t();
//       Rprintf("--bk(%d): ",bk.Storage()); //printCVector(bk); 
      Fk=Fsample.SubMatrix(n,n,(k-1)*m*(m*p+1)+1,k*m*(m*p+1)).t(); pFk=Fk.Store(); 
//       Rprintf("--Fk(%d): ",Fk.Storage()); //printCVector(Fk); 

      for(i=1;i<=m;i++) xik(i)=sqrt(xisample(n,(k-1)*m+i)); 
//       Rprintf("--xik(%d)/sqrtxik(%d) defined\n",m,m); 

      // Compute A0/A0^-1/sqrtwish for regime k
      A0=b2a(bk,Ui); 
      //Rprintf("--A0(%d,%d):",m,m); //printMatrix(A0); 
      A0i=A0.i(); 
      //Rprintf("--A0^-1(%d,%d):",m,m); //printMatrix(A0i); 
      sqrtwish=(A0*xik).i(); 
      //Rprintf("--sqrtwish(%d,%d):",m,m); //printMatrix(sqrtwish); 

      // Compute beta vector 
      fmat.ReSize(m,m*p+1); fmat<<pFk; fmat=fmat.t(); 
      fmat=(fmat.Rows(1,m*p)*A0i).t(); bvec=fmat.AsColumn(); 
//       Rprintf("--fmat(%d,%d):",m,m*p+1); printMatrix(fmat); 
//       Rprintf("bvec_%d:", n); printCVector(bvec);
      
      // Compute IRF 
      impulse.Row(n)=irf_var_from_beta(sqrtwish.t(), bvec, ns).t(); 
      if (!(n%N210pct))
	Rprintf("Regime %d: Monte Carlo IRF %d percent complete (Iteration %d)\n",k,++pctct*10,n);
    }

    // Create and class Robj for impulses, load into IRFlist
    PROTECT(IRFtmp=C2R3D(impulse,IRFdims)); 
    setclass(IRFtmp,"mc.irf.BSVAR"); SET_VECTOR_ELT(IRFlist, k-1, IRFtmp); 
    UNPROTECT(1); 
  }
  UNPROTECT(5); 
  return IRFlist;
}
