#include "QRD.h"

// Empty constructor 
QRD::QRD(){_n=-1;}

// 1) Perform QR Decomp via LINPACK 'dqrdc2' subroutine 
// 2) Populate QRD member data with decomposed QR matrix and
//    info about decomp routine 
QRD::QRD(const Matrix& A, const double tol=1E-12) 
{
  // Initialize member data and allocate heap memory
  _n=A.Nrows(); _p=A.Ncols(); _tol=tol; 
  _qr=C2F(A);
  _pivot=(int*)R_alloc(_p,sizeof(int)); 
  _tau=(double*)R_alloc(_p,sizeof(double)); 
  _work=(double*)R_alloc(_p*2,sizeof(double));
  for(int i=0;i<_p;i++) 
  {
	_pivot[i]=i+1; 
	_tau[i]=0; 
	_work[i]=0; 
	_work[i+_p]=0; 
  }

  // LINPACK QR factorization via householder transformations
  F77_CALL(dqrdc2)(_qr, &_n, &_n, &_p, &_tol, &_rank, _tau, _pivot, _work);
}

// Return Generalized Orthogonal Factor Q(n,n)
ReturnMatrix QRD::Q() const
{
  // Pull member data, allocate memory and initialize Y/Q
  int i, j, ct=0, n=_n, rank=_rank; double *y, *q;
  y=(double *)R_alloc(n*n,sizeof(double)); 
  q=(double *)R_alloc(n*n,sizeof(double)); 
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
    { 
	y[ct++] = (i==j) ? 1 : 0; 
	q[ct-1]=y[ct-1]; 
    }

  // Compute full orthogonal factor Q(n,n) 
  F77_CALL(dqrqy)(_qr, &n, &rank, _tau, y, &n, q); 
  
  // Convert FORTRAN to C and return
  Matrix Qout=F2C(q,n,n);
  return Qout.ForReturn(); 
}

// Return Generalized Upper Triangular Factor R(n,p)
ReturnMatrix QRD::R() const
{
  int i, j, minnp=_n<_p?_n:_p;
  Matrix R(minnp,_p); 
  R=0.0; 
  double *r=R.Store();
  for(i=0;i<minnp;i++) 
    for(j=i;j<_p;j++) 
	r[i*_p+j] = _qr[j*_n+i];
  return R.ForReturn(); 
}

// QR Solve Routine via LINPACK subroutine dqrcf
ReturnMatrix QRD::Solve(const Matrix &Y) const
{
  // Catch problems with QRD/Y 
  if(_n<1 || _rank!=_p || Y.Nrows()!=_n)
  { 
	Matrix t(1,1); 
	t=-1; 
	return t.ForReturn(); 
  }

  // Assign input values and allocate memory
  int i, n=_n, k=_rank, ny=Y.Ncols(), info[1]; 
  info[0]=1;
  double *y=C2F(Y), *coef=(double*)R_alloc(k*ny,sizeof(double)); 
  for(i=0;i<k*ny;i++) 
	coef[i]=0.0; 

  // Call FORTRAN routine dqrcf
  F77_CALL(dqrcf)(_qr, &n, &k, _tau, y, &ny, coef, info);
  
  // Check result
  if(info[0]!=0){
    Rprintf("\nLINPACK subroutine 'dqrcf' found exact singularity\n");
    Matrix t(1,1); 
    t=-1;  
    return t.ForReturn(); 
  }
  
  // Build return object, clean up and return
  Matrix out=F2C(coef,k,ny); 
  return out.ForReturn(); 
} 

// Accessor functions for decomposition information

// (packed) QR Matrix
ReturnMatrix QRD::QR() const { return F2C(_qr,_n,_p); }

// column pivoting info
int QRD::pivot(int idx) const { return _pivot[idx]; }
int* QRD::pivot() const { return _pivot; }

// qr$aux: householder reflection multiplier 
double QRD::tau(int idx) const { return _tau[idx]; }
double* QRD::tau() const { return _tau; } 

int QRD::rank() const { return _rank; }     // rank
int QRD::nrow() const { return _n; }        // nrow(QR)
int QRD::ncol() const { return _p; }        // ncol(QR)
