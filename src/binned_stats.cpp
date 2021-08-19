#include "binned_values.h"

#include <stdexcept>
#include <unordered_map>
#include <map>
#include <vector>

#include <cpp11.hpp>
#include <Rmath.h>

#include <boost/functional/hash.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/extended_p_square.hpp>

#include "logging.h"

typedef std::vector<int> var_set_t;
typedef int r_index_t;
typedef std::pair<r_index_t, r_index_t> vargroup_id_t;
typedef std::unordered_map<vargroup_id_t, var_set_t, boost::hash<vargroup_id_t>> vargroup_map_t;

using namespace cpp11;

namespace bacc = boost::accumulators;

//??? Probability that random X-distributed variable would
//??? be less or equal than the given value.
//??? Gaussian kernel smoothing is used to represent the
//??? distribution of X variable.
//???
//??? @X samples of X random variable
//??? @y value to compare with
//??? @nsteps the number of segments to divide the range of X value into
//??? @bandwidth the gaussian smoothing kernel bandwidth, defaults to the segment size, 0 disables smoothing
//??? @return P(X<=y)
double ProbabilityLessSmoothed(
    doubles X,
    double y,
    int   nsteps = 100,
    double bandwidth = na<double>()
){
    if ( X.size() == 0 ) {
        throw std::length_error( "X is empty" );
    }
    ImportedValues xvals( X );

    if ( ( xvals.val_max - xvals.val_min ) > nsteps * std::numeric_limits<double>::epsilon() ) {
        return BinnedValues( xvals, xvals.defaultBinWidth( nsteps ), true ).probabilityLessOrEqual( y, bandwidth );
    } else {
        // degenerated
        LOG_DEBUG1( "Degenerated distribution" );
        if ( is_na( bandwidth ) ) {
            if ( fabs( xvals.val_max - y ) <= nsteps * std::numeric_limits<double>::epsilon() ) {
                LOG_DEBUG1( "Near zero" );
                return ( 0.5 );
            } else {
                return ( xvals.val_max <= y ? 1.0 : 0.0 );
            }
        } else {
            return ( Rf_pnorm5( y, xvals.val_max, bandwidth, 1, 0 ) );
        }
    }
}

//??? Probability that random X-distributed variable would
//??? be less or equal than the Y-distributed variable.
//??? Gaussian kernel smoothing is used to represent the
//??? distribution of X-Y variable.
//???
//??? @X samples of X random variable
//??? @Y samples of Y random variable
//??? @nsteps the number of segments to divide the range of X-Y value into
//??? @bandwidth the gaussian smoothing kernel bandwidth, defaults to the segment size, 0 disables smoothing
//??? @return P(X<Y)
[[cpp11::register]]
double ProbabilityLessSmoothed(
    doubles  X,
    doubles  Y,
    int   nsteps,
    double bandwidth
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

    return ( diffBins.probabilityLessOrEqual( 0.0, bandwidth ) );
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

//??? Average probability that random X-distributed variable would
//??? be less or equal than the Y-distributed variable.
//??? Gaussian kernel smoothing is used to represent the
//??? distribution of X-Y variable.
//??? The probability is averaged across the set of X and Y experiments
//???
//??? @X matrix of samples of X random variable, columns are different experiments, rows are MCMC iterations
//??? @Y matrix of samples of Y random variable, columns are different experiments, rows are MCMC iterations
//??? @Deltas 3-D array of adjustments to add to Y when comparing i-th experiment of X and j-th experiment of Y
//??? @nsteps the number of segments to divide the range of X-Y value into
//??? @maxBandwidth the maximum gaussian smoothing kernel bandwidth, 0 disables smoothing
//??? @matchIterations if true, X[i]-Y[i] is used, otherwise X[i]-X[j], where i and j are independent iterations
//??? @return E(X-Y), D(X-Y), P(X<Y)
[[cpp11::register]]
writable::data_frame DifferenceStatistics(
    doubles_matrix<by_column>  X,
    doubles_matrix<by_column>  Y,
    doubles  Deltas,
    int   nsteps,
    double maxBandwidth,
    bool  matchIterations
){
    if ( X.ncol() == 0 || X.nrow() == 0 ) {
        throw std::length_error( "X matrix is empty" );
    }
    if ( Y.ncol() == 0 || Y.nrow() == 0 ) {
        throw std::length_error( "Y matrix is empty" );
    }
    if ( !is_na( maxBandwidth ) && (maxBandwidth < 0.0) ) {
        throw std::invalid_argument( "maxBandwidth cannot be negative" );
    }

    integers deltaDims(Deltas.attr("dim"));
    if ( deltaDims.size() != 3 ) {
        throw std::invalid_argument( "Deltas is not 3D array" );
    }
    if ( deltaDims[2] != X.nrow() ) {
        throw std::length_error( "Number of rows in X and 3rd dimension in Deltas does not match" );
    }
    if ( deltaDims[2] != Y.nrow() ) {
        throw std::length_error( "Number of rows in Y and 3rd dimension in Deltas does not match" );
    }
    list_of<strings> deltaDimNames(list(Deltas.attr("dimnames")));
    LOG_DEBUG2( "Loading delta 1st dimension names" );
    name_index_map_t deltaXNames = nameToIndexMap( deltaDimNames[0] );
    LOG_DEBUG2( "Loading delta 2nd dimension names" );
    name_index_map_t deltaYNames = nameToIndexMap( deltaDimNames[1] );
    std::size_t deltaStride = deltaDims[0] * deltaDims[1]; // stride of vectors in the 3rd dimension
    std::size_t ndraws = deltaDims[2];

    LOG_DEBUG2( "Loading X column names" );
    strings xColNames(list_of<strings>(list(X.attr("dimnames"))[1]));
    LOG_DEBUG2( "Loading Y column names" );
    strings yColNames(list_of<strings>(list(Y.attr("dimnames"))[1]));

    std::vector<double> xvals_merged;
    std::vector<double> yvals_adj_merged;
    xvals_merged.reserve( X.ncol() * Y.ncol() * ndraws );
    yvals_adj_merged.reserve( X.ncol() * Y.ncol() * ndraws );

    for ( int i = 0; i < X.ncol(); i++ ) {
        ImportedValues xvals( X[i] );

        name_index_map_t::const_iterator iColIt = deltaXNames.find( xColNames[i] );
        if ( iColIt == deltaXNames.end() ) {
            THROW_EXCEPTION(std::domain_error,
                            "Column '%s' not found in Deltas 1st dimension", as_cpp<const char *>(xColNames[i]));
        }

        for ( int j = Y.ncol()-1; j >= 0; j-- ) {
            LOG_DEBUG2("Processing distribution of %s vs %s",
                       as_cpp<const char *>(xColNames[i]), as_cpp<const char *>(yColNames[j]));
            name_index_map_t::const_iterator jColIt = deltaYNames.find( yColNames[j] );
            if ( jColIt == deltaYNames.end() ) {
                THROW_EXCEPTION(std::domain_error,
                                "Column '%s' not found in Deltas 2nd dimension", as_cpp<const char *>(yColNames[j]));
            }
            xvals_merged.insert( xvals_merged.end(), xvals.values.cbegin(), xvals.values.cend() );

            auto yCol(Y[j]);
            std::size_t deltaOffset = iColIt->second + jColIt->second * deltaDims[0];
            for ( std::size_t k = 0; k < ndraws; k++ ) {
                yvals_adj_merged.push_back( yCol[k] + Deltas[ deltaOffset + k * deltaStride ] );
            }
        }
    }
    LOG_DEBUG1( "Calculating probability that difference is non-positive" );
    double prob;
    double diff_mean;
    double diff_var;
    double bw;
    if ( matchIterations ) {
        // calculate per-iteration deltas
        std::vector<double> diffs( yvals_adj_merged.size() );
        for ( std::size_t k = 0; k < diffs.size(); k++ ) {
            diffs[k] = xvals_merged[k] - yvals_adj_merged[k];
        }
        ImportedValues diff_vals( diffs );
        BinnedValues diffBins( diff_vals, diff_vals.defaultBinWidth( nsteps ) );
        bw = maxBandwidth == 0.0 ? 0.0 : diffBins.bw_nrd();
        if ( !is_na( maxBandwidth ) && (bw > maxBandwidth) ) {
            bw = maxBandwidth;
        }
        prob = diffBins.probabilityLessOrEqual( 0.0, bw );
        diff_mean = diffBins.average();
        diff_var = diffBins.variance();
    } else {
        // bandwidth is calculated from the global difference of X values vs Y values
        BinnedValues globalDiffBins = BinnedValues::difference( ImportedValues( xvals_merged ),
                                                                ImportedValues( yvals_adj_merged ), nsteps );
        bw = maxBandwidth == 0.0 ? 0.0 : globalDiffBins.bw_nrd();
        if ( !is_na( maxBandwidth ) && (bw > maxBandwidth) ) {
            bw = maxBandwidth;
        }
        // calculate probability for every X column and all Y columns
        // NOTE it's not correct to calculate the probability using globalDiffBins,
        // because every X-Y column pair has its own adjustment
        std::size_t chunk_size = Y.ncol() * ndraws;
        std::vector<double> probs;
        double max_prob = 0;
        double sum_mean = 0;
        double sum_var = 0;
        for ( int i = 0; i < X.ncol(); i++ ) {
            BinnedValues diffBins = BinnedValues::difference( ImportedValues( std::vector<double>( xvals_merged.begin() + i * chunk_size,
                                                                              xvals_merged.begin() + (i+1) * chunk_size ) ),
                                                              ImportedValues( std::vector<double>( yvals_adj_merged.begin() + i * chunk_size,
                                                                              yvals_adj_merged.begin() + (i+1) * chunk_size ) ), nsteps );
            probs.push_back( diffBins.probabilityLessOrEqual( 0.0, bw ) );
            if ( probs.size() == 1 || max_prob < probs.back() ) {
                max_prob = probs.back();
            }
            sum_mean += diffBins.average();
            sum_var += diffBins.variance();
        }
        LOG_DEBUG( "Calculating average probability" );
        if ( max_prob == 0.0 ) {
            // replace 0 P-value with minimal double, so that log-transform is still possible
            prob = std::numeric_limits<double>::min();
        } else {
            double sum_prob = 0.0;
            for ( std::size_t i = 0; i < probs.size(); i++ ) {
                sum_prob += probs[i] / max_prob;
            }
            prob = max_prob * sum_prob / probs.size();
        }
        diff_mean = sum_mean / X.ncol();
        diff_var = sum_var / X.ncol();
    }
    LOG_DEBUG("P(X <= Y) = %g", prob);
    return ( writable::data_frame{
        "mean"_nm = diff_mean,
        "sd"_nm = sqrt( diff_var ),
        "prob_less"_nm = prob,
        "bw"_nm = bw});
}

#if 0 // unsupported default params
writable::data_frame ContrastStatistics(
        doubles   X,
        integers  col2group_col,
        integers  col2group_group,
        doubles_matrix<by_row>  contrastXvargroup,
        doubles  contrast_offsets,
        int   nsteps = 100,
        double maxBandwidth = na<double>(),
        doubles  quant_probs = writable::doubles{0.025, 0.25, 0.50, 0.75, 0.975}
#endif

//??? Calculates average probabilities that given contrasts would be
//??? be less or equal than zero.
//??? Gaussian kernel smoothing is used to represent the contrast distribution.
//??? The probability is averaged across all possible combinations of
//??? the columns of the relevant groups.
//???
//??? @var2group_var index of the variable (X column) in many-to-many var <-> group map
//??? @var2group_group index of the group in many-to-many var <-> group map
//??? @var2group_contrast index of the contrast in many-to-many var <-> group map (group contents could be contrast-specific)
//??? @contrastXvargroup contrast matrix, rows are contrasts, columns are variable groups
//??? @contrast_offsets vector of contrast offsets, i.e. the reported difference is not X-Y, it's X-Y+offset
//??? @X matrix of samples of X random variable, columns are different variables, rows are MCMC draws
//??? @nsteps the number of segments to divide the range of X values into
//??? @maxBandwidth the maximum Gaussian smoothing kernel bandwidth, 0 disables smoothing
//??? @return E(X), D(X), P(X>=0), P(X<=0)
[[cpp11::register]]
data_frame ContrastStatistics(
    doubles   X,
    integers  var2group_var,
    integers  var2group_group,
    integers  var2group_contrast,
    doubles_matrix<by_row>  contrastXvargroup,
    doubles  contrast_offsets,
    int   nsteps,
    double maxBandwidth,
    function summaryfun
){
    R_xlen_t ncontrasts = contrastXvargroup.nrow();
    R_xlen_t ngroups = contrastXvargroup.ncol();
    LOG_DEBUG1("ncontrasts=%d ngroups=%d", ncontrasts, ngroups);
    if ( contrastXvargroup.ncol() == 0 || contrastXvargroup.nrow() == 0 ) {
        THROW_EXCEPTION( std::length_error, "contrastXvargroup matrix is empty" );
    }
    if ( X.size() == 0 ) {
        THROW_EXCEPTION( std::length_error, "X draws are empty" );
    }
    if ( var2group_var.size() != var2group_group.size() ) {
        THROW_EXCEPTION(std::length_error, "var2group vars length (%d) doesn't match the groups length (%d)",
                        var2group_var.size(), var2group_group.size());
    }
    if ( var2group_var.size() != var2group_contrast.size() ) {
        THROW_EXCEPTION(std::length_error, "var2group vars length (%d) doesn't match the contrasts length (%d)",
                        var2group_var.size(), var2group_contrast.size());
    }
    if ( !is_na( maxBandwidth ) && (maxBandwidth < 0.0) ) {
        THROW_EXCEPTION(std::domain_error, "maxBandwidth cannot be negative");
    }
    if ( contrast_offsets.size() != ncontrasts ) {
        THROW_EXCEPTION(std::length_error, "contrast_offsets size (%d) doesn't match the rows of contrastXvargroup (%d)",
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
    LOG_DEBUG1("nmcmc_iters=%d nchains=%d nvars=%d", niters, nchains, nvars);

    LOG_DEBUG1("Initializing vargroup map");
    vargroup_map_t vargroups;
    for ( std::size_t i = 0; i < var2group_group.size(); i++ ) {
        r_index_t group_ix = var2group_group[i];
        if ( group_ix == 0 ) {
            THROW_EXCEPTION(std::domain_error, "Group index should be in 1:%d range, %d found", ngroups, group_ix);
        } else {
            --group_ix;
        }
        r_index_t contrast_ix = var2group_contrast[i];
        if ( contrast_ix == 0 ) {
            THROW_EXCEPTION(std::domain_error, "Contrast index should be in 1:%d range, %d found", ncontrasts, contrast_ix);
        } else {
            --contrast_ix;
        }
        r_index_t var_ix = var2group_var[i];
        if ( var_ix == 0 || var_ix > nvars ) {
            THROW_EXCEPTION(std::domain_error, "Variable index should be in 1:%d range, %d found", nvars, var_ix);
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
    LOG_DEBUG1( "Calculating the probabilities that contrast is non-positive and non-negative" );
    writable::doubles contrast_draws{na<double>()}; // init to something as otherwise it is null
    writable::integers contrast_draws_dims{(int)niters, (int)nchains, 1};
    contrast_draws.attr("class") = writable::strings{"array"};
    for ( int contr_ix = 0; contr_ix < ncontrasts; ++contr_ix ) {
        // calculate how many conditions are associated with current contrast
        LOG_DEBUG1("Calculating contrast #%d", contr_ix);
        std::vector<vargroup_map_t::const_iterator> contr_groups;
        std::size_t nperm = 1;
        auto group_weights = contrastXvargroup[contr_ix];

        bool empty_groups = false;
        for ( R_xlen_t group_ix = 0; group_ix < ngroups; ++group_ix ) {
            const double w = group_weights[group_ix];
            if (w != 0.0) {
                LOG_DEBUG2("contrastXvargroup[%d, %d]=%g", contr_ix, group_ix, w);
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
        LOG_DEBUG1( "  " << nperm << " permutations" );
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
