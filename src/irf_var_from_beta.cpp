#include "MSBVARcpp.h"

//extern "C"
ReturnMatrix irf_var_from_beta(const Matrix& A0, const ColumnVector& bvec, const int nsteps)
{
  int m=A0.Ncols();
  int p=bvec.Storage()/(m*m), i,j;

  Matrix bmat(m,m*p), AR(m,m*p); bmat = 0.0; AR=0.0; bmat << bvec.Store();
  for(i=1;i<=p;i++)
    AR.SubMatrix(1,m,1+m*(i-1),m*i) = bmat.SubMatrix(1,m,1+m*(i-1),m*i).t();

  Matrix mhat=irf_var_mhat(AR, nsteps, A0), tmp(nsteps,m);
  ColumnVector IRF;
  for(i=1;i<=m;i++){
    tmp=0.0;
    for(j=1;j<=nsteps;j++) tmp.Row(j) = mhat.Column(i+(j-1)*m).t();
    if(i==1) IRF = (tmp.t()).AsColumn(); else IRF &= (tmp.t()).AsColumn();
  }
  IRF.Release(); return IRF.ForReturn();
}

//extern "C"
ReturnMatrix irf_var_mhat(const Matrix& AR, const int nsteps, const Matrix& A0)
{
  int i, j, m=AR.Nrows(), p=AR.Ncols()/m;
  int maxp = (nsteps > p) ? nsteps : p, minp = (nsteps < p) ? nsteps : p;
  Matrix Bmat(m,m*maxp); Bmat=0.0;
  Bmat.SubMatrix(1,m,1,m*minp) = AR.SubMatrix(1,m,1,m*minp);

  // IRF Matrix and Identification condition
  Matrix Mhat(m,m*nsteps); Mhat = 0.0; Mhat.SubMatrix(1,m,1,m)=A0;

  // Compute IRF
  for(i=2; i<=nsteps; i++)
    for(j=1; j<=i-1; j++)
      Mhat.SubMatrix(1,m,1+m*(i-1),m*i) =
	Mhat.SubMatrix(1,m,1+m*(i-1),m*i) +
	(Mhat.SubMatrix(1,m,1+m*(i-j-1),m*(i-j)) *
	 Bmat.SubMatrix(1,m,1+m*(j-1),m*j));

  Mhat.Release(); return Mhat.ForReturn();
}
