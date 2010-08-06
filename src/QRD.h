#ifndef QRD_H
#define QRD_H 

#include "MSBVARcpp.h"

// QRD:  Class interface for QR-related decompositions,
//       factorizations, and least squares solutions.

class QRD{ 
  int _n;        // nrow(A) 
  int _p;        // ncol(A) 
  int _rank;     // rank
  int *_pivot;   // pivot information  
  double _tol;   // tolerance  
  double *_tau;  // householder reflection info  
  double *_work; // work vector for fortran routines 
  double *_qr;   // QR decomposition in packed form 
  //  double *_x;    // Column major input matrix 
  
 public:
  QRD(); //~QRD();            

  // A=QR decomp via LINPACK subroutine dqrdc2
  QRD(const Matrix&, const double); 
  
  // QR utility functions
  ReturnMatrix QR() const;     // Get packed QR(n,p)
  ReturnMatrix Q() const;      // Unpack factor Q(n,n)
  ReturnMatrix R() const;      // Unpack factor R(n,p)
  ReturnMatrix Solve(const Matrix &) const;  // Ax=b 
  
  // Private data member accessor functions
  int* pivot() const;    // _pivot array
  int pivot(int) const;  // _pivot[idx]

  double* tau() const;      // _tau array (qr$aux)
  double tau(int) const;    // _tau[idx]

  int rank() const;      // _rank
  int nrow() const;      // _n
  int ncol() const;      // _p
};

#endif

// Body file: QRD.cpp 
