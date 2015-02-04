/*
*  wrapper for calling R's random number generator from
*  the original FORTRAN code for FFBS.f.  Passes the unif value from R
*  to the Fortran function MVUNI which is then going to bingen() for
*  the backward sampler.
*
*  Modeled after the same in mvtnorm pkg
*
*/
#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>

void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }
double F77_SUB(unifrnd)(void) { return unif_rand(); }

