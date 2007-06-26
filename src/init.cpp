#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "MSBVARcpp.h"

extern "C" {

//   SEXP irf_var(SEXP, SEXP, SEXP, SEXP);
//   SEXP mc_irf_var(SEXP, SEXP, SEXP);
//   SEXP msbsvar_irf(SEXP, SEXP, SEXP); 
//   SEXP drawA0(SEXP, SEXP, SEXP, SEXP, SEXP);
//   SEXP gibbsA0(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//   SEXP log_marg_A0k(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
//   SEXP mc_irf_bsvar(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);  

//   SEXP a2b_test(SEXP, SEXP, SEXP, SEXP); 
//   SEXP irf_var_from_beta_debug(SEXP, SEXP, SEXP);
//   SEXP BHLK_FILTER(SEXP, SEXP, SEXP);
//   SEXP BHLK_SMOOTHER(SEXP, SEXP);
//   SEXP SSencode(SEXP);
//   SEXP SSdecode(SEXP);

  // Define .Call methods 
  static R_CallMethodDef callMethods[] = {
    {"irf.var.cpp", (DL_FUNC) &irf_var, 4},
    {"mc.irf.var.cpp", (DL_FUNC) &mc_irf_var, 3},
    {"mc.irf.bsvar.cpp", (DL_FUNC) &mc_irf_bsvar, 12}, 
    {"irf.msbsvar.cpp", (DL_FUNC) &msbsvar_irf, 3},
    {"drawA0.cpp", (DL_FUNC) &drawA0, 5},
    {"gibbsA0.cpp",(DL_FUNC) &gibbsA0, 6},
    {"log.marginal.A0k",(DL_FUNC) &log_marg_A0k, 9},
    {NULL, NULL, 0}
  };

//     {"a2b.cpp",(DL_FUNC) &a2b_test, 4}, 
//     {"irf.var.from.beta.debug", (DL_FUNC) &irf_var_from_beta_debug, 3},
//     {"BHLK.filter.cpp", (DL_FUNC) &BHLK_FILTER, 3},
//     {"BHLK.smoother.cpp", (DL_FUNC) &BHLK_SMOOTHER, 2},
//     {"SS.encode", (DL_FUNC) &SSencode, 1},
//     {"SS.decode", (DL_FUNC) &SSdecode, 1},
  
  void R_init_MSBVAR(DllInfo *info)
  {  
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);  
  }
}
