#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Lapack.h>  // LAPACK entry points
#include <R_ext/Linpack.h> // LINPACK entry points
#include <R_ext/Applic.h>
#include <math.h>
//#include <stdbool.h>
//#include <stdint.h>        // uint_32t class type
#include <limits.h>        // POSIX standardization

/* #include <string> */
/* #include <vector> */
using namespace std;

#include "newmatap.h"         // matrix apps + newmat.h
#include "newmatrc.h"         // matrix row/column functions
using namespace NEWMAT;

//extern "C" {

#define LOG2PI log(2*M_PI)

#include "A0_W.h"          // Class for storing gibbs draws
#include "QRD.h"           // Class for QR routines
//#include "SS.h"            // State Space Class
//#include "newSS.h"	   // New State Space Class

  // R object creation/manipulation functions
  void setdims(SEXP, int, int*);
  int *getdims(SEXP);
  void setclass(SEXP, const char*);
  SEXP makeList(SEXP *, char**);
  SEXP listElt(SEXP, const char*);

  // R to C object translation
  ReturnMatrix R2CRV(SEXP);
  ReturnMatrix R2CCV(SEXP);
  ReturnMatrix R2Cmat(SEXP, int, int);
  ReturnMatrix R2Cmat2(SEXP);
  ReturnMatrix R2Carr(SEXP);

  // C to R object translation
  SEXP C2Rmat(const Matrix&);
  SEXP C2R3D(const Matrix&, int *);
  SEXP C2R3D2(const Matrix&, int *);
  SEXP C2Rint(int*);
  SEXP C2Rdouble(double *);
  SEXP C2Rdoublemat(double *, int, int);

  // C <==> FORTRAN object translation
  ReturnMatrix F2C(double *, int, int);
  double *C2F(const Matrix &);

  // Fast, in-place row-major/col-major translation
/*   void rm2cm_double(double *, int, int); */
/*   void cm2rm_double(double *, int, int); */
/*   void rm2cm_Matrix(Matrix &mat);  */

  // Simple math util functions for class Matrix
  ReturnMatrix absmat(const Matrix&);
  ReturnMatrix sqrtVec(const Matrix&);
  ReturnMatrix cumprod(const Matrix&);
  ReturnMatrix colsums(const Matrix&);

  // Extra Matrix Subset Functions
  ReturnMatrix diag(const Matrix&);
  ReturnMatrix rep(const ColumnVector&, int);

  // Matrix Print Function Overloads
  void printMatrix(const Matrix&);
  void printRVector(const RowVector&);
  void printCVector(const ColumnVector&);

  // Simple math util functions for numerical stability
  Real sqrtHyp(const Real&, const Real&);

  // Functions to draw from different distributions
  ReturnMatrix rnorms(int);
  ReturnMatrix rnorms_mat(int, int);
  ReturnMatrix rwish(const Matrix&, int);
  ReturnMatrix dmvnorm(Matrix&, Matrix&);

  // C versions of MSBVAR hidden functions
  ReturnMatrix irf_var_from_beta(const Matrix&, const ColumnVector&, const int);
  ReturnMatrix irf_var_mhat(const Matrix& , const int, const Matrix&);
  ReturnMatrix b2a(const ColumnVector&, SEXP);
  ReturnMatrix a2b(const Matrix&, SEXP, int*, int*);
  ReturnMatrix norm_svar(const Matrix&, const Matrix&, const int, int*);
  ReturnMatrix drawA0cpp(const Matrix&, const class UTobj&, const int, const int*, class Wobj&);

  // Posterior fit helper function for log marginal A0k computation
  SEXP log_marg_A0k(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

  // C versions of MSBVAR functions declared in the namespace
  SEXP irf_var(SEXP, SEXP, SEXP, SEXP);
  SEXP mc_irf_var(SEXP, SEXP, SEXP);
  SEXP mc_irf_bsvar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP msbsvar_irf(SEXP, SEXP, SEXP);
  SEXP drawA0(SEXP, SEXP, SEXP, SEXP, SEXP);
  SEXP gibbsA0(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

  //} // end extern "C"
