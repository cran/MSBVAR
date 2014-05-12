#include "MSBVARcpp.h"

//extern "C"
SEXP irf_var(SEXP ar_coefs, SEXP ar_dims, SEXP nsteps, SEXP A0)
{
  SEXP B, mhat, names, output;
  int m=INTEGER(ar_dims)[0], p=INTEGER(ar_dims)[2],ns=INTEGER(nsteps)[0];
  int i, j, k, maxp = (ns > p) ? ns : p, minp = (ns < p) ? ns : p;

  // Define Matrix for AR coefs
  // Fill Matrix tmp with m*m elements of R obj
  // Assign transpose of tmp to appropriate submatrix to fix
  // row-major/col-major difference between R/C++

  Real *rar = REAL(ar_coefs);
  Matrix Bmat(m,m*maxp); Bmat = 0.0;
  Matrix Btmp(m,m);
  R_len_t ct=0;
  for(i=1;i<=minp;i++){
    Btmp = 0.0;
    for(j=1;j<=m;j++)
      for(k=1;k<=m;k++)
	Btmp(j,k) = rar[ct++];
    Bmat.SubMatrix(1,m,1+m*(i-1),m*i) = Btmp.t();
  }

  // Array to hold IRF
  Matrix Mhat(m,m*ns); Mhat = 0.0;

  // Identification condition
  Mhat.SubMatrix(1,m,1,m) = R2Cmat(A0,m,m);

  // Compute IRF
  for(i=2; i<=ns; i++)
    for(j=1; j<=i-1; j++)
      Mhat.SubMatrix(1,m,1+m*(i-1),m*i) =
	Mhat.SubMatrix(1,m,1+m*(i-1),m*i) +
	(Mhat.SubMatrix(1,m,1+m*(i-j-1),m*(i-j)) *
	 Bmat.SubMatrix(1,m,1+m*(j-1),m*j));

  // Create R objects 'B' (AR coefs) and 'mhat' (IRF)
  int dims[]={m,m,ns};
  PROTECT(B = C2R3D(Bmat, dims));
  PROTECT(mhat = C2R3D(Mhat, dims));

  // Create R list object 'output'
  PROTECT(output = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(output, 0, B);
  SET_VECTOR_ELT(output, 1, mhat);

  // Set output names
  PROTECT(names = allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("B"));
  SET_STRING_ELT(names, 1, mkChar("mhat"));
  setAttrib(output, R_NamesSymbol, names);
  UNPROTECT(1);

  // Set output class
  setclass(output,"irf.VAR");

  UNPROTECT(3);
  return output;
} // end irf_var()

