#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "MSBVARcpp.h"

extern "C" {

  // Define .Call methods 
  static R_CallMethodDef callMethods[] = {
    {"irf.var.cpp", (DL_FUNC) &irf_var, 4},
    {"mc.irf.var.cpp", (DL_FUNC) &mc_irf_var, 3},
    {"mc.irf.bsvar.cpp", (DL_FUNC) &mc_irf_bsvar, 12}, 
    {"irf.msbsvar.cpp", (DL_FUNC) &msbsvar_irf, 3},
    {"drawA0.cpp", (DL_FUNC) &drawA0, 5},
    {"gibbsA0.cpp",(DL_FUNC) &gibbsA0, 6},
    {"log.marginal.A0k",(DL_FUNC) &log_marg_A0k, 9},
    //    {"irf.var.from.beta", (DL_FUNC) &irf_var_from_beta, 3},
    {NULL, NULL, 0}
  };
  
  void R_init_MSBVAR(DllInfo *info)
  {  
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);  
  }
}
