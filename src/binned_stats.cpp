#include "binned_values.h"

#include <stdexcept>
#include <unordered_map>
#include <map>
#include <vector>
#include <tuple>
#include <numeric>

#include <cpp11.hpp>
#include <Rmath.h>

#include <boost/functional/hash.hpp>

#include "logging.h"
#include "pvalue_utils.h"

typedef std::vector<int> var_set_t;
typedef int r_index_t;
typedef std::pair<r_index_t, r_index_t> vargroup_id_t;
typedef std::unordered_map<vargroup_id_t, var_set_t, boost::hash<vargroup_id_t>> vargroup_map_t;

using namespace cpp11;

typedef std::map<r_string, std::size_t> name_index_map_t;

name_index_map_t nameToIndexMap( const strings& strs )
{
    name_index_map_t res;
    for ( std::size_t i = 0; i < (std::size_t)strs.size(); i++ ) {
        res.insert( std::make_pair( strs[i], i ) );
    }
    return ( res );
}

#if 0 // unsupported default params
writable::data_frame ContrastStatistics(
        doubles   X,
        integers  col2group_col,
        integers  col2group_group,
        doubles_matrix<by_row>  vargroupXcontrast,
        doubles  contrast_offsets,
        int   nsteps = 100,
        double maxBandwidth = na<double>(),
        doubles  quant_probs = writable::doubles{0.025, 0.25, 0.50, 0.75, 0.975},
        double mlog10pvalue_threshold = 10.0,
        double mlog10pvalue_hard_threshold_factor = 3.0
#endif

// Calculates average probabilities that given contrasts would be
// be less or equal than zero using *contrast_calculator* function.
// The probability is averaged across all possible combinations of
// the columns of the relevant groups.
//
// @param var2group_var index of the variable (X column) in many-to-many var <-> group map
// @param var2group_group index of the group in many-to-many var <-> group map
// @param var2group_contrast index of the contrast in many-to-many var <-> group map (group contents could be contrast-specific)
// @param vargroupXcontrast contrast matrix, rows are contrasts, columns are variable groups
// @param contrast_offsets vector of contrast offsets, i.e. the reported difference is not X-Y, it's X-Y+offset
// @param mlog10pvalue_threshold if `-log10(pvalue)` is above this threshold, compresses p-value using `pvalue_sqrt_compress()`
//
// @return E(X), D(X), P(X>=0), P(X<=0)
//
// @seealso [ContrastStatistics_draws()]
template<typename T>
data_frame ContrastStatistics(
    T         contrast_calculator,
    integers  var2group_var,
    integers  var2group_group,
    integers  var2group_contrast,
    doubles_matrix<by_column>  vargroupXcontrast,
    doubles   contrast_offsets,
    double    mlog10pvalue_threshold,
    double    mlog10pvalue_hard_threshold_factor
){
  R_xlen_t nvars = contrast_calculator.nvars();
  R_xlen_t ncontrasts = vargroupXcontrast.ncol();
  R_xlen_t ngroups = vargroupXcontrast.nrow();
  LOG_DEBUG1("ncontrasts=%d ngroups=%d", ncontrasts, ngroups);
  if (ncontrasts == 0 || ngroups == 0) {
    THROW_EXCEPTION(std::length_error, "vargroupXcontrast matrix is empty");
  }
  if (var2group_var.size() != var2group_group.size()) {
    THROW_EXCEPTION(std::length_error, "var2group vars length (%ld) doesn't match the groups length (%ld)",
                    var2group_var.size(), var2group_group.size());
  }
  if (var2group_var.size() != var2group_contrast.size()) {
    THROW_EXCEPTION(std::length_error, "var2group vars length (%ld) doesn't match the contrasts length (%ld)",
                    var2group_var.size(), var2group_contrast.size());
  }
  if (contrast_offsets.size() != ncontrasts) {
    THROW_EXCEPTION(std::length_error, "contrast_offsets size (%ld) doesn't match the rows of vargroupXcontrast (%ld)",
                    contrast_offsets.size(), ncontrasts);
  }

  LOG_DEBUG1("Initializing vargroup map");
  vargroup_map_t vargroups;
  for (R_xlen_t i = 0; i < var2group_group.size(); i++) {
    r_index_t group_ix = var2group_group[i];
    if (group_ix == 0) {
      THROW_EXCEPTION(std::domain_error, "Group index should be in 1:%ld range, %d found", ngroups, group_ix);
    } else {
      --group_ix;
    }
    r_index_t contrast_ix = var2group_contrast[i];
    if (contrast_ix == 0) {
      THROW_EXCEPTION(std::domain_error, "Contrast index should be in 1:%ld range, %d found", ncontrasts, contrast_ix);
    } else {
      --contrast_ix;
    }
    r_index_t var_ix = var2group_var[i];
    if (var_ix == 0 || var_ix > nvars) {
      THROW_EXCEPTION(std::domain_error, "Variable index should be in 1:%ld range, %d found", nvars, var_ix);
    } else {
      --var_ix;
    }
    const auto cXg = std::make_pair(contrast_ix, group_ix);
    vargroup_map_t::iterator groupIt = vargroups.find(cXg);
    if (groupIt == vargroups.end()) {
      vargroups.insert(groupIt, std::make_pair(cXg, var_set_t(1, var_ix)));
    } else {
      groupIt->second.push_back(var_ix);
    }
  }

  writable::integers index_contrast;
  writable::doubles probs_nonpos;
  writable::doubles probs_nonneg;
  writable::doubles nperms;

  writable::list extra_dfs{};

  // generate contrast samples for all relevant columns combinations
  LOG_DEBUG1("Calculating the probabilities that contrast is non-positive and non-negative");
  for ( int contr_ix = 0; contr_ix < ncontrasts; ++contr_ix ) {
    // calculate how many conditions are associated with current contrast
    LOG_DEBUG1("Calculating contrast #%d", contr_ix);
    std::vector<vargroup_map_t::const_iterator> contr_groups;
    std::size_t nperm = 1;
    auto group_weights = vargroupXcontrast[contr_ix];

    bool empty_groups = false;
    for ( R_xlen_t group_ix = 0; group_ix < ngroups; ++group_ix ) {
      const double w = group_weights[group_ix];
      if (w != 0.0) {
        LOG_DEBUG2("vargroupXcontrast[%d, %d]=%g", contr_ix, group_ix, w);
        const auto cXg = vargroup_id_t(contr_ix, group_ix);
        vargroup_map_t::const_iterator groupIt = vargroups.find(cXg);
        if ( groupIt != vargroups.end() ) {
          contr_groups.push_back(groupIt);
          nperm *= groupIt->second.size();
        } else {
          empty_groups = true;
          break;
          //THROW_EXCEPTION(std::invalid_argument, "No variables for group #" << (group_ix+1) << " found");
        }
      }
    }
    if (empty_groups) {
      LOG_DEBUG1(" skipping contrast #%d computation since some groups were empty", contr_ix);
      continue;
    }
    LOG_DEBUG1("  %d permutations", nperm);
    if (nperm == 0) { // should not happen as checked before
      THROW_EXCEPTION(std::invalid_argument, "No variable combinations for contrast #%d", contr_ix);
    }
    const double contr_offset = contrast_offsets[contr_ix];
    std::vector<var_set_t::const_iterator> var_its;
    std::transform(contr_groups.begin(), contr_groups.end(), std::back_inserter(var_its),
                   [](const vargroup_map_t::const_iterator& it){ return it->second.cbegin(); });
    contrast_calculator.init_contrast(contr_offset, nperm);
    for (std::size_t perm_i = 0; perm_i < nperm; ++perm_i) {
      LOG_DEBUG2("Processing contrast permutation #%d of %d", perm_i, nperm);
      // generate contrast samples for the next permutation of current contrast
      contrast_calculator.begin_permutation();
      for (std::size_t group_i = 0; group_i < contr_groups.size(); ++group_i) {
        const auto groupIt = contr_groups[group_i];
        const auto group_ix = groupIt->first.second;
        const double w = group_weights[group_ix];
        const r_index_t var_ix = *(var_its[group_i]);
        LOG_DEBUG2("Processing var #%d of group #%d", var_ix, group_ix+1);
        contrast_calculator.add_variable(var_ix, w);
      }
      // increment vars-in-groups enumeration
      for (std::size_t group_i = 0; group_i < contr_groups.size(); ++group_i) {
        auto&& it = var_its[group_i];
        const auto& vars = contr_groups[group_i]->second;
        it++;
        if (it == vars.cend()) {
          it = vars.cbegin(); // start over
        } else {
          break; // inc done
        }
      }
      if ((LOG_LEVEL >= 3) && (var_its.size() > 1)) {
        std::ostringstream permvars;
        for (const auto var_it : var_its) {
          permvars << " " << (*var_it);
        }
        LOG_DEBUG("perm #%d: %s", perm_i, permvars.str().c_str());
      }
      contrast_calculator.end_permutation();
    }

    // calculate the contrast
    double prob_nonpos, prob_nonneg;
    sexp extra_df;
    std::tie(prob_nonpos, prob_nonneg, extra_df) = contrast_calculator.calculate();
    index_contrast.push_back(contr_ix + 1);
    probs_nonpos.push_back(prob_nonpos);
    probs_nonneg.push_back(prob_nonneg);
    nperms.push_back(nperm);
    extra_dfs.push_back(extra_df);
    LOG_DEBUG2("%d-th contrast done", contr_ix);
  }
  LOG_DEBUG2("contrasts done");

  if (!is_na(mlog10pvalue_threshold)) {
    LOG_DEBUG2("compressing p-values");
    for (auto it = probs_nonneg.begin(); it != probs_nonneg.end(); ++it) {
      *it = pvalue_sqrt_compress(*it, mlog10pvalue_threshold,
                                 mlog10pvalue_hard_threshold_factor);
    }
    for (auto it = probs_nonpos.begin(); it != probs_nonpos.end(); ++it) {
      *it = pvalue_sqrt_compress(*it, mlog10pvalue_threshold,
                                 mlog10pvalue_hard_threshold_factor);
    }
  }

  LOG_DEBUG2("composing resulting data.frame");
  static auto bindcols = package("dplyr")["bind_cols"];
  static auto bindrows = package("dplyr")["bind_rows"];
  data_frame res_df = writable::data_frame{
      "__contrast_ix__"_nm = index_contrast,
      "prob_nonneg"_nm = probs_nonneg,
      "prob_nonpos"_nm = probs_nonpos,
      "nperms"_nm = nperms
    };
  R_xlen_t nextra_dfs = std::accumulate(extra_dfs.cbegin(), extra_dfs.cend(), 0,
                                    [](R_xlen_t n, sexp df) -> R_xlen_t
                                    { return !Rf_isNull(df) ? n+1 : n; });
  LOG_DEBUG2("%ld contrast(s) have extra data frames", nextra_dfs);
  if (nextra_dfs > 0 && nextra_dfs < res_df.nrow()) {
    cpp11::warning("Some summary data frames are null, summary is ignored");
  }
  return nextra_dfs == res_df.nrow() ?
    as_cpp<data_frame>(bindcols(res_df, bindrows(extra_dfs))) :
    res_df;
}

// contrast_calculator based on MCMC draws
class DrawsContrastCalculator {
private:
  doubles  draws_;
  sexp     summaryfun_;
  const int      nsteps_;
  const double   maxBandwidth_;

  R_xlen_t niters_;
  R_xlen_t nchains_;
  R_xlen_t nvars_;

  double offset_;
  size_t nperms_;

  size_t perm_i_;
  writable::doubles contrast_draws_;
  writable::integers contrast_draws_dims_;

public:
  DrawsContrastCalculator(doubles draws, sexp summaryfun,
                          int nsteps, double maxBandwidth) :
    draws_(draws), summaryfun_(summaryfun),
    nsteps_(nsteps), maxBandwidth_(maxBandwidth),
    niters_(0), nchains_(0), nvars_(0),
    offset_(na<double>()), nperms_(0),
    perm_i_(0),
    contrast_draws_{na<double>()}, // init to something as otherwise it is null
    contrast_draws_dims_{0, 0, 1}
  {
    if (draws.size() == 0) {
      THROW_EXCEPTION(std::length_error, "draws are empty");
    }
    if (!is_na(maxBandwidth) && (maxBandwidth < 0.0)) {
      THROW_EXCEPTION(std::domain_error, "maxBandwidth cannot be negative");
    }
    integers draw_dims(draws.attr("dim"));
    if (draw_dims.size() != 3) {
      THROW_EXCEPTION(std::length_error, "draws is not a 3D (iterations x chains x variables) array");
    }
    niters_ = draw_dims[0];
    nchains_ = draw_dims[1];
    nvars_ = draw_dims[2];
    LOG_DEBUG1("nmcmc_iters=%ld nchains=%ld nvars=%ld", niters_, nchains_, nvars_);
    contrast_draws_.attr("class") = writable::strings{"array"};
    contrast_draws_dims_[1] = nchains_;
  }

  R_xlen_t nvars() const {
    return nvars_;
  }

  void init_contrast(double offset, size_t nperms) {
    offset_ = offset;
    nperms_ = nperms;

    // reset draws
    contrast_draws_.attr("dim") = R_NilValue; // temporary unset before resizing
    contrast_draws_.resize(niters_ * nchains_ * nperms_);
    contrast_draws_dims_[0] = niters_ * nperms_;
    contrast_draws_.attr("dim") = contrast_draws_dims_; // set back the updated dims
    std::fill(contrast_draws_.begin(), contrast_draws_.end(), 0.0);
    perm_i_ = 0;
  }

  void begin_permutation() {
    if (is_na<double>(offset_)) throw std::runtime_error("contrast is not initialized");
    if (perm_i_ >= nperms_) THROW_EXCEPTION(std::domain_error, "permutation #%ld out of range (%ld)",
                                            perm_i_+1, nperms_);
  }
  void end_permutation() {
    if (is_na<double>(offset_)) throw std::runtime_error("contrast is not initialized");
    perm_i_++;
    if (perm_i_ > nperms_) THROW_EXCEPTION(std::domain_error, "permutation #%ld out of range (%ld)",
                                           perm_i_+1, nperms_);
  }
  void add_variable(r_index_t var_ix, double w) {
    if (is_na<double>(offset_)) throw std::runtime_error("contrast is not initialized");

    const auto perm_end = contrast_draws_.begin() + perm_i_ + contrast_draws_.size();
    auto srcit = draws_.cbegin() + var_ix * niters_ * nchains_;
    for (auto dstit = contrast_draws_.begin() + perm_i_;
         dstit != perm_end; ++srcit, (dstit += nperms_)) {
      (*dstit) += w * (*srcit);
    }
  }
  std::tuple<double, double, data_frame> calculate() const {
    if (is_na<double>(offset_)) throw std::runtime_error("contrast is not initialized");
    if (perm_i_ != nperms_) throw std::runtime_error("incomplete permutations");

    // calculate the contrast
    ImportedValues contrast_vals(contrast_draws_);
    BinnedValues contrast_bins(contrast_vals, contrast_vals.defaultBinWidth(nsteps_));

    double cur_bw = maxBandwidth_ == 0.0 ? 0.0 : contrast_bins.bw_nrd();
    if (!is_na(maxBandwidth_) && (cur_bw > maxBandwidth_)) {
      cur_bw = maxBandwidth_;
    }
    double prob_nonpos = contrast_bins.probabilityLessOrEqual(offset_, cur_bw);
    double prob_nonneg = contrast_bins.probabilityGreaterOrEqual(offset_, cur_bw);

    //LOG_RCOUT("step=" << contrast_bins.step << " min=" << contrast_bins.val_min << " max=" << contrast_bins.val_max);
    //Rcpp::Rcout << "bins=";
    //for ( int i = 0; i < contrast_bins.bins.size(); i++ ) {
    //    Rcpp::Rcout << contrast_bins.bins[i];
    //    Rcpp::Rcout << " ,";
    //}
    //Rcpp::Rcout << "\n";
    //LOG_RCOUT("mean_direct=" << std::accumulate(all_contrast_samples.begin(), all_contrast_samples.end(), 0.0) / all_contrast_samples.size());
    //LOG_RCOUT("sd_direct=" << sqrt(Rcpp::sugar::Var(all_contrast_samples)));

    static auto bindcols = package("dplyr")["bind_cols"];

    const data_frame bw_df = writable::data_frame{
                                "bandwidth"_nm = cur_bw,
                                "bin_width"_nm = contrast_bins.step};
    const data_frame summary_df = !Rf_isNull(summaryfun_) ?
        as_cpp<data_frame>(bindcols(
            as_cpp<data_frame>((as_cpp<function>(summaryfun_))(contrast_draws_)),
            bw_df)) :
        bw_df;

    return std::make_tuple(prob_nonpos, prob_nonneg, summary_df);
  }
};

// Calculates average probabilities that given contrasts would be
// be less or equal than zero.
// Gaussian kernel smoothing is used to represent the contrast distribution.
// The probability is averaged across all possible combinations of
// the columns of the relevant groups.
//
// @param draws matrix of samples of random variables, columns are different variables, rows are MCMC draws
// @param var2group_var index of the variable (X column) in many-to-many var <-> group map
// @param var2group_group index of the group in many-to-many var <-> group map
// @param var2group_contrast index of the contrast in many-to-many var <-> group map (group contents could be contrast-specific)
// @param vargroupXcontrast contrast matrix, rows are contrasts, columns are variable groups
// @param contrast_offsets vector of contrast offsets, i.e. the reported difference is not X-Y, it's X-Y+offset
// @param nsteps the number of segments to divide the range of X values into
// @param maxBandwidth the maximum Gaussian smoothing kernel bandwidth, 0 disables smoothing
// @param mlog10pvalue_threshold if `-log10(pvalue)` is above this threshold, compresses p-value using `pvalue_sqrt_compress()`
//
// @return E(X), D(X), P(X>=0), P(X<=0)
[[cpp11::register]]
data_frame ContrastStatistics_draws(
    doubles   draws,
    integers  var2group_var,
    integers  var2group_group,
    integers  var2group_contrast,
    doubles_matrix<by_column>  vargroupXcontrast,
    doubles   contrast_offsets,
    int       nsteps,
    double    maxBandwidth,
    sexp      summaryfun,
    double    mlog10pvalue_threshold,
    double    mlog10pvalue_hard_threshold_factor
){
    return ContrastStatistics(DrawsContrastCalculator(draws, summaryfun, nsteps, maxBandwidth),
                              var2group_var, var2group_group, var2group_contrast,
                              vargroupXcontrast, contrast_offsets,
                              mlog10pvalue_threshold, mlog10pvalue_hard_threshold_factor);
}
