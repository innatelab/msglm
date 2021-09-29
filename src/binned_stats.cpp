#include "binned_values.h"

#include <stdexcept>
#include <unordered_map>
#include <map>
#include <vector>

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

// Probability that random X-distributed variable would
// be less or equal than the given value.
// Gaussian kernel smoothing is used to represent the
// distribution of X variable.
//
// @param X samples of X random variable
// @param y value to compare with
// @param nsteps the number of segments to divide the range of X value into
// @param bandwidth the gaussian smoothing kernel bandwidth, defaults to the segment size, 0 disables smoothing
// @param mlog10_threshold if `-log10(prob)` is above this threshold, compresses *prob* using `pvalue_sqrt_compress()`
// @return P(X<=y)
double ProbabilityLessSmoothed(
    doubles X,
    double y,
    int   nsteps = 100,
    double bandwidth = na<double>(),
    double mlog10_threshold = 10.0,
    double mlog10_hard_threshold_factor = 3.0
){
    if ( X.size() == 0 ) {
        throw std::length_error( "X is empty" );
    }
    ImportedValues xvals( X );
    double res = na<double>();

    if ( ( xvals.val_max - xvals.val_min ) > nsteps * std::numeric_limits<double>::epsilon() ) {
        res = BinnedValues( xvals, xvals.defaultBinWidth( nsteps ), true ).probabilityLessOrEqual( y, bandwidth );
    } else {
        // degenerated
        LOG_DEBUG1( "Degenerated distribution" );
        if ( is_na( bandwidth ) ) {
            if ( fabs( xvals.val_max - y ) <= nsteps * std::numeric_limits<double>::epsilon() ) {
                LOG_DEBUG1( "Near zero" );
                res = 0.5;
            } else {
                res = xvals.val_max <= y ? 1.0 : 0.0;
            }
        } else {
            res = Rf_pnorm5( y, xvals.val_max, bandwidth, 1, 0 );
        }
    }
    if (!is_na(mlog10_threshold)) {
        res = pvalue_sqrt_compress(res, mlog10_threshold, mlog10_hard_threshold_factor);
    }
    return res;
}

// Probability that random X-distributed variable would
// be less or equal than the Y-distributed variable.
// Gaussian kernel smoothing is used to represent the
// distribution of X-Y variable.
//
// @param X samples of X random variable
// @param Y samples of Y random variable
// @param nsteps the number of segments to divide the range of X-Y value into
// @param bandwidth the gaussian smoothing kernel bandwidth, defaults to the segment size, 0 disables smoothing
// @param mlog10_threshold if `-log10(prob)` is above this threshold, compresses *prob* using `pvalue_sqrt_compress()`
// @return P(X<Y)
[[cpp11::register]]
double ProbabilityLessSmoothed(
    doubles  X,
    doubles  Y,
    int   nsteps,
    double bandwidth,
    double mlog10_threshold = 10.0,
    double mlog10_hard_threshold_factor = 3.0
){
    if ( X.size() == 0 ) {
        throw std::length_error( "X is empty" );
    }
    if ( Y.size() == 0 ) {
        throw std::length_error( "Y is empty" );
    } else if ( Y.size() == 1 ) { // compare X samples against single Y value
        return ProbabilityLessSmoothed(X, Y[0], nsteps, bandwidth);
    }
    ImportedValues xvals( X );
    ImportedValues yvals( Y );
    BinnedValues diffBins = BinnedValues::difference( xvals, yvals, nsteps );

    double res = diffBins.probabilityLessOrEqual(0.0, bandwidth);
    if (!is_na(mlog10_threshold)) {
        res = pvalue_sqrt_compress(res, mlog10_threshold, mlog10_hard_threshold_factor);
    }
    return res;
}

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
// be less or equal than zero.
// Gaussian kernel smoothing is used to represent the contrast distribution.
// The probability is averaged across all possible combinations of
// the columns of the relevant groups.
//
// @param var2group_var index of the variable (X column) in many-to-many var <-> group map
// @param var2group_group index of the group in many-to-many var <-> group map
// @param var2group_contrast index of the contrast in many-to-many var <-> group map (group contents could be contrast-specific)
// @param vargroupXcontrast contrast matrix, rows are contrasts, columns are variable groups
// @param contrast_offsets vector of contrast offsets, i.e. the reported difference is not X-Y, it's X-Y+offset
// @param X matrix of samples of X random variable, columns are different variables, rows are MCMC draws
// @param nsteps the number of segments to divide the range of X values into
// @param maxBandwidth the maximum Gaussian smoothing kernel bandwidth, 0 disables smoothing
// @param mlog10pvalue_threshold if `-log10(pvalue)` is above this threshold, compresses p-value using `pvalue_sqrt_compress()`
// @return E(X), D(X), P(X>=0), P(X<=0)
[[cpp11::register]]
data_frame ContrastStatistics(
    doubles   X,
    integers  var2group_var,
    integers  var2group_group,
    integers  var2group_contrast,
    doubles_matrix<by_column>  vargroupXcontrast,
    doubles  contrast_offsets,
    int   nsteps,
    double maxBandwidth,
    function summaryfun,
    double mlog10pvalue_threshold,
    double mlog10pvalue_hard_threshold_factor
){
    R_xlen_t ncontrasts = vargroupXcontrast.ncol();
    R_xlen_t ngroups = vargroupXcontrast.nrow();
    LOG_DEBUG1("ncontrasts=%d ngroups=%d", ncontrasts, ngroups);
    if ( ncontrasts == 0 || ngroups == 0 ) {
        THROW_EXCEPTION( std::length_error, "vargroupXcontrast matrix is empty" );
    }
    if ( X.size() == 0 ) {
        THROW_EXCEPTION( std::length_error, "X draws are empty" );
    }
    if ( var2group_var.size() != var2group_group.size() ) {
        THROW_EXCEPTION(std::length_error, "var2group vars length (%ld) doesn't match the groups length (%ld)",
                        var2group_var.size(), var2group_group.size());
    }
    if ( var2group_var.size() != var2group_contrast.size() ) {
        THROW_EXCEPTION(std::length_error, "var2group vars length (%ld) doesn't match the contrasts length (%ld)",
                        var2group_var.size(), var2group_contrast.size());
    }
    if ( !is_na( maxBandwidth ) && (maxBandwidth < 0.0) ) {
        THROW_EXCEPTION(std::domain_error, "maxBandwidth cannot be negative");
    }
    if ( contrast_offsets.size() != ncontrasts ) {
        THROW_EXCEPTION(std::length_error, "contrast_offsets size (%ld) doesn't match the rows of vargroupXcontrast (%ld)",
                        contrast_offsets.size(), ncontrasts);
    }

    integers Xdims(X.attr("dim"));
    if ( Xdims.size() != 3 ) {
        THROW_EXCEPTION(std::length_error, "X is not a 3D (draws x chains x variables) array" );
    }
    R_xlen_t niters = Xdims[0];
    R_xlen_t nchains = Xdims[1];
    R_xlen_t nvars = Xdims[2];
    R_xlen_t ndraws = niters * nchains;
    LOG_DEBUG1("nmcmc_iters=%ld nchains=%ld nvars=%ld", niters, nchains, nvars);

    LOG_DEBUG1("Initializing vargroup map");
    vargroup_map_t vargroups;
    for ( std::size_t i = 0; i < var2group_group.size(); i++ ) {
        r_index_t group_ix = var2group_group[i];
        if ( group_ix == 0 ) {
            THROW_EXCEPTION(std::domain_error, "Group index should be in 1:%ld range, %d found", ngroups, group_ix);
        } else {
            --group_ix;
        }
        r_index_t contrast_ix = var2group_contrast[i];
        if ( contrast_ix == 0 ) {
            THROW_EXCEPTION(std::domain_error, "Contrast index should be in 1:%ld range, %d found", ncontrasts, contrast_ix);
        } else {
            --contrast_ix;
        }
        r_index_t var_ix = var2group_var[i];
        if ( var_ix == 0 || var_ix > nvars ) {
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
    writable::doubles prob_nonpos;
    writable::doubles prob_nonneg;
    writable::doubles bw;
    writable::doubles bin_width;
    writable::doubles nperms;
    writable::doubles niterations;

    sexp postsummary_df;
    auto rbind = package("base")["rbind"];

    // generate contrast samples for all relevant columns combinations
    LOG_DEBUG1("Calculating the probabilities that contrast is non-positive and non-negative");
    writable::doubles contrast_draws{na<double>()}; // init to something as otherwise it is null
    writable::integers contrast_draws_dims{(int)niters, (int)nchains, 1};
    contrast_draws.attr("class") = writable::strings{"array"};
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
        std::vector<var_set_t::const_iterator> var_its;
        std::transform(contr_groups.begin(), contr_groups.end(), std::back_inserter(var_its),
                       [](const vargroup_map_t::const_iterator& it){ return it->second.cbegin(); });
        contrast_draws.attr("dim") = R_NilValue; // temporary unset before resizing
        contrast_draws.resize(ndraws * nperm);
        contrast_draws_dims[0] = niters * nperm;
        contrast_draws.attr("dim") = contrast_draws_dims; // set back the updated dims
        std::fill(contrast_draws.begin(), contrast_draws.end(), 0.0);
        for (std::size_t perm_i = 0; perm_i < nperm; ++perm_i) {
            LOG_DEBUG2("Processing contrast permutation #%d of %d", perm_i, nperm);
            // generate contrast samples for the next permutation of current contrast
            const auto perm_begin = contrast_draws.begin() + perm_i;
            const auto perm_end = contrast_draws.begin() + perm_i + ndraws * nperm;
            for (std::size_t group_i = 0; group_i < contr_groups.size(); ++group_i) {
                const auto groupIt = contr_groups[group_i];
                const auto group_ix = groupIt->first.second;
                const double w = group_weights[group_ix];
                const r_index_t var_ix = *(var_its[group_i]);
                auto srcit = X.cbegin() + var_ix * ndraws;
                LOG_DEBUG2("Processing var #%d of group #%d", var_ix, group_ix+1);
                for (auto dstit = perm_begin;
                     dstit != perm_end; ++srcit, (dstit += nperm)) {
                    (*dstit) += w * (*srcit);
                }
            }
            // increment vars-in-groups enumeration
            for ( std::size_t group_i = 0; group_i < contr_groups.size(); ++group_i ) {
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
        }

        // calculate the contrast
        ImportedValues contrast_vals( contrast_draws );
        BinnedValues contrast_bins( contrast_vals, contrast_vals.defaultBinWidth( nsteps ) );

        double contr_offset = contrast_offsets[contr_ix];
        double cur_bw = maxBandwidth == 0.0 ? 0.0 : contrast_bins.bw_nrd();
        if ( !is_na( maxBandwidth ) && (cur_bw > maxBandwidth) ) {
            cur_bw = maxBandwidth;
        }
        index_contrast.push_back(contr_ix + 1);
        prob_nonpos.push_back(contrast_bins.probabilityLessOrEqual(contr_offset, cur_bw));
        prob_nonneg.push_back(contrast_bins.probabilityGreaterOrEqual(contr_offset, cur_bw));
        bw.push_back(cur_bw);
        bin_width.push_back(contrast_bins.step);
        nperms.push_back(nperm);
        niterations.push_back(niters);

        const auto cur_summary_df = summaryfun(contrast_draws);
        postsummary_df = SEXP(postsummary_df) != R_NilValue ? rbind(postsummary_df, cur_summary_df) : cur_summary_df;
        //LOG_RCOUT("step=" << contrast_bins.step << " min=" << contrast_bins.val_min << " max=" << contrast_bins.val_max);
        //Rcpp::Rcout << "bins=";
        //for ( int i = 0; i < contrast_bins.bins.size(); i++ ) {
        //    Rcpp::Rcout << contrast_bins.bins[i];
        //    Rcpp::Rcout << " ,";
        //}
        //Rcpp::Rcout << "\n";
        //LOG_RCOUT("mean_direct=" << std::accumulate(all_contrast_samples.begin(), all_contrast_samples.end(), 0.0) / all_contrast_samples.size());
        //LOG_RCOUT("sd_direct=" << sqrt(Rcpp::sugar::Var(all_contrast_samples)));
        //LOG_DEBUG( "P[>=0]=" << prob_nonneg <<
        //           " P[<=0]=" << prob_nonpos );
        LOG_DEBUG2("%d-th contrast done", contr_ix);
    }
    LOG_DEBUG2("contrasts done");

    if (!is_na(mlog10pvalue_threshold)) {
        LOG_DEBUG2("compressing p-values");
        for (auto it = prob_nonneg.begin(); it != prob_nonneg.end(); ++it) {
            *it = pvalue_sqrt_compress(*it, mlog10pvalue_threshold, mlog10pvalue_hard_threshold_factor);
        }
        for (auto it = prob_nonpos.begin(); it != prob_nonpos.end(); ++it) {
            *it = pvalue_sqrt_compress(*it, mlog10pvalue_threshold, mlog10pvalue_hard_threshold_factor);
        }
    }
    LOG_DEBUG2("composing resulting data.frame");
    return as_cpp<data_frame>(package("base")["cbind"](writable::data_frame{
           "__contrast_ix__"_nm = index_contrast
        },
        postsummary_df,
        writable::data_frame{
            "prob_nonneg"_nm = prob_nonneg,
            "prob_nonpos"_nm = prob_nonpos,
            "bin_width"_nm = bin_width,
            "bw"_nm = bw,
            "nperms"_nm = nperms,
            "niters"_nm = niterations
        }));
}
