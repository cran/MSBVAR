#include "MSBVARcpp.h"

// ================================================================
// A0 normalization rules
// ================================================================
//
// method=0:  DistanceMLA
//
//   Multiply the inverse of A0 by the mode of A0. If the sign of
//   a diagonal element in the result is negative, flip the sign of
//   the elements of the corresponding column of A0.
//
//   Increase switch count.
//
// method=1:  DistanceMLAhat
//
//   Premultiply A0 by the inverse of the mode of A0.If the sign of
//   a diagonal element in the result is negative, flip the sign of
//   the elements of the corresponding column of A0.
//
//   Increase switch count.
//
// method=2:  Euclidean
//
//   Find the column sums of the sqaure of the difference between A0
//   and its mode and -A0 and its mode. Flip the sign of the columns
//   of A0 where the column sums for -A0 are smaller.
//
//   Increase switch count.
//
// method=3:	PositiveDiagA
//
//   Flip the sign of the columns of A0 where the diagonal element
//   of that column is negative.
//
// method=4:  PositiveDiagAinv
//
//   Flip the sign of the columns of A0 where the diagnoal element
//   of the corresponding column of the inverse of A0 is negative.
//
//   Increase switch count.
//
// ================================================================
// ================================================================

//extern "C"
ReturnMatrix norm_svar(const Matrix& A0in,
				  const Matrix& A0mode,
				  const int norm_method,
				  int *switchct)
{
  int i, n=A0in.Ncols(); Matrix A0=A0in, tmp, tmpn; RowVector csums, csumsn;

  switch (norm_method){
  case 0:
    tmp = A0.i()*A0mode;
    for(i=1;i<=n;i++)	if(tmp(i,i)<0){A0.Column(i)=-A0.Column(i); switchct++;}
    break;

  case 1:
    tmp = A0mode.i()*A0;
    for(i=1;i<=n;i++) if(tmp(i,i)<0){A0.Column(i)=-A0.Column(i); switchct++;}
    break;

  case 2:
    tmp=A0-A0mode; tmp*=tmp; tmpn=-A0-A0mode; tmpn*=tmpn;
    csums=colsums(tmp); csumsn=colsums(tmpn);
    for(i=1;i<=n;i++)
	if(csumsn(i)<csums(i)){A0.Column(i)=-A0.Column(i); switchct++;}
    break;

  case 3:
    for(i=1;i<=n;i++)	if(A0(i,i)<0){A0.Column(i)=-A0.Column(i); switchct++;}
    break;

  case 4:
    tmp=A0.i();
    for(i=1;i<=n;i++)	if(tmp(i,i)<0){A0.Column(i)=-A0.Column(i); switchct++;}
    break;
  }

  A0.Release(); return A0.ForReturn();
}
