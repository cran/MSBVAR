#include "MSBVARcpp.h"

//extern "C"
ReturnMatrix rnorms(int num)
{
  ColumnVector norms(num);
  for(int i=1; i<=num; i++) norms(i) = rnorm(0,1);
  norms.Release(); return norms.ForReturn();
}

//extern "C"
ReturnMatrix rnorms_mat(int nr, int nc)
{
  int i, j; Matrix norms(nr,nc);
  for(i=1;i<=nr;i++) for(j=1;j<=nc;j++) norms(i,j)=rnorm(0,1);
  norms.Release(); return norms.ForReturn();
}

//extern "C"

ReturnMatrix rwish(const Matrix& Sigma, int df)
{
  int p = Sigma.Nrows();
  SymmetricMatrix tmp; tmp << Sigma;
  Matrix ZS = rnorms(p*df).AsMatrix(p,df).t() * Cholesky(tmp).t();
  Matrix wish = ZS.t()*ZS;
  wish.Release(); return wish.ForReturn();
}
