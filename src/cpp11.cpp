// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// binned_stats.cpp
double ProbabilityLessSmoothed(doubles  X, doubles  Y, int   nsteps, double bandwidth, double mlog10_threshold, double mlog10_hard_threshold_factor);
extern "C" SEXP _msglm_ProbabilityLessSmoothed(SEXP X, SEXP Y, SEXP nsteps, SEXP bandwidth, SEXP mlog10_threshold, SEXP mlog10_hard_threshold_factor) {
  BEGIN_CPP11
    return cpp11::as_sexp(ProbabilityLessSmoothed(cpp11::as_cpp<cpp11::decay_t<doubles >>(X), cpp11::as_cpp<cpp11::decay_t<doubles >>(Y), cpp11::as_cpp<cpp11::decay_t<int  >>(nsteps), cpp11::as_cpp<cpp11::decay_t<double>>(bandwidth), cpp11::as_cpp<cpp11::decay_t<double>>(mlog10_threshold), cpp11::as_cpp<cpp11::decay_t<double>>(mlog10_hard_threshold_factor)));
  END_CPP11
}
// binned_stats.cpp
writable::data_frame DifferenceStatistics(doubles_matrix<by_column>  X, doubles_matrix<by_column>  Y, doubles  Deltas, int   nsteps, double maxBandwidth, bool  matchIterations, double mlog10pvalue_threshold, double mlog10pvalue_hard_threshold_factor);
extern "C" SEXP _msglm_DifferenceStatistics(SEXP X, SEXP Y, SEXP Deltas, SEXP nsteps, SEXP maxBandwidth, SEXP matchIterations, SEXP mlog10pvalue_threshold, SEXP mlog10pvalue_hard_threshold_factor) {
  BEGIN_CPP11
    return cpp11::as_sexp(DifferenceStatistics(cpp11::as_cpp<cpp11::decay_t<doubles_matrix<by_column> >>(X), cpp11::as_cpp<cpp11::decay_t<doubles_matrix<by_column> >>(Y), cpp11::as_cpp<cpp11::decay_t<doubles >>(Deltas), cpp11::as_cpp<cpp11::decay_t<int  >>(nsteps), cpp11::as_cpp<cpp11::decay_t<double>>(maxBandwidth), cpp11::as_cpp<cpp11::decay_t<bool >>(matchIterations), cpp11::as_cpp<cpp11::decay_t<double>>(mlog10pvalue_threshold), cpp11::as_cpp<cpp11::decay_t<double>>(mlog10pvalue_hard_threshold_factor)));
  END_CPP11
}
// binned_stats.cpp
data_frame ContrastStatistics(doubles   X, integers  var2group_var, integers  var2group_group, integers  var2group_contrast, doubles_matrix<by_column>  vargroupXcontrast, doubles  contrast_offsets, int   nsteps, double maxBandwidth, function summaryfun, double mlog10pvalue_threshold, double mlog10pvalue_hard_threshold_factor);
extern "C" SEXP _msglm_ContrastStatistics(SEXP X, SEXP var2group_var, SEXP var2group_group, SEXP var2group_contrast, SEXP vargroupXcontrast, SEXP contrast_offsets, SEXP nsteps, SEXP maxBandwidth, SEXP summaryfun, SEXP mlog10pvalue_threshold, SEXP mlog10pvalue_hard_threshold_factor) {
  BEGIN_CPP11
    return cpp11::as_sexp(ContrastStatistics(cpp11::as_cpp<cpp11::decay_t<doubles  >>(X), cpp11::as_cpp<cpp11::decay_t<integers >>(var2group_var), cpp11::as_cpp<cpp11::decay_t<integers >>(var2group_group), cpp11::as_cpp<cpp11::decay_t<integers >>(var2group_contrast), cpp11::as_cpp<cpp11::decay_t<doubles_matrix<by_column> >>(vargroupXcontrast), cpp11::as_cpp<cpp11::decay_t<doubles >>(contrast_offsets), cpp11::as_cpp<cpp11::decay_t<int  >>(nsteps), cpp11::as_cpp<cpp11::decay_t<double>>(maxBandwidth), cpp11::as_cpp<cpp11::decay_t<function>>(summaryfun), cpp11::as_cpp<cpp11::decay_t<double>>(mlog10pvalue_threshold), cpp11::as_cpp<cpp11::decay_t<double>>(mlog10pvalue_hard_threshold_factor)));
  END_CPP11
}
// cloud_stats.cpp
writable::data_frame CloudsOverlap(doubles_matrix<by_column>  X, doubles_matrix<by_column>  Y, double maxBandwidth, int   nsteps, doubles quant_probs);
extern "C" SEXP _msglm_CloudsOverlap(SEXP X, SEXP Y, SEXP maxBandwidth, SEXP nsteps, SEXP quant_probs) {
  BEGIN_CPP11
    return cpp11::as_sexp(CloudsOverlap(cpp11::as_cpp<cpp11::decay_t<doubles_matrix<by_column> >>(X), cpp11::as_cpp<cpp11::decay_t<doubles_matrix<by_column> >>(Y), cpp11::as_cpp<cpp11::decay_t<double>>(maxBandwidth), cpp11::as_cpp<cpp11::decay_t<int  >>(nsteps), cpp11::as_cpp<cpp11::decay_t<doubles>>(quant_probs)));
  END_CPP11
}
// pvalue_combine.cpp
double GammaMethodPValuesCombine(doubles pvalues, double stt);
extern "C" SEXP _msglm_GammaMethodPValuesCombine(SEXP pvalues, SEXP stt) {
  BEGIN_CPP11
    return cpp11::as_sexp(GammaMethodPValuesCombine(cpp11::as_cpp<cpp11::decay_t<doubles>>(pvalues), cpp11::as_cpp<cpp11::decay_t<double>>(stt)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_msglm_CloudsOverlap",             (DL_FUNC) &_msglm_CloudsOverlap,              5},
    {"_msglm_ContrastStatistics",        (DL_FUNC) &_msglm_ContrastStatistics,        11},
    {"_msglm_DifferenceStatistics",      (DL_FUNC) &_msglm_DifferenceStatistics,       8},
    {"_msglm_GammaMethodPValuesCombine", (DL_FUNC) &_msglm_GammaMethodPValuesCombine,  2},
    {"_msglm_ProbabilityLessSmoothed",   (DL_FUNC) &_msglm_ProbabilityLessSmoothed,    6},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_msglm(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}