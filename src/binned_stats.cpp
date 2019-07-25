// [[Rcpp::depends(BH)]]

#include "binned_values.h"

#include <Rcpp.h>

typedef std::vector<int> experiment_set_t;
#if __has_include(<unordered_map>)
    #include <unordered_map>
    typedef std::unordered_map<std::size_t, experiment_set_t> experiment_map_t;
#elif __has_include(<tr1/unordered_map>)
    #include <tr1/unordered_map>
    typedef std::tr1::unordered_map<std::size_t, experiment_set_t> experiment_map_t;
#endif

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/extended_p_square.hpp>

#include "logging.h"

namespace bacc = boost::accumulators;

//??? Probability that random X-distributed variable would
//??? be less or equal than zero.
//??? Gaussian kernel smoothing is used to represent the
//??? distribution of X-Y variable.
//???
//??? @X samples of X random variable
//??? @nsteps the number of segments to divide the range of X value into
//??? @bandwidth the gaussian smoothing kernel bandwidth, defaults to the segment size, 0 disables smoothing
//??? @return P(X<=0)
// [[Rcpp::export]]
double ProbabilityLessZeroSmoothed(
    const Rcpp::NumericVector& X,
    int   nsteps = 100,
    double bandwidth = NA_REAL
){
    if ( X.size() == 0 ) {
        throw Rcpp::exception( "X is empty" );
    }
    ImportedValues xvals( X );
    if ( ( xvals.val_max - xvals.val_min ) > nsteps * std::numeric_limits<double>::epsilon() ) {
        return BinnedValues( xvals, xvals.defaultBinWidth( nsteps ), true ).probabilityNonPositive( bandwidth );
    } else {
        // degenerated
        LOG_DEBUG1( "Degenerated distribution" );
        if ( R_IsNA( bandwidth ) ) {
            if ( fabs( xvals.val_max ) <= nsteps * std::numeric_limits<double>::epsilon() ) {
                LOG_DEBUG1( "Near zero" );
                return ( 0.5 );
            } else {
                return ( xvals.val_max <= 0.0 ? 1.0 : 0.0 );
            }
        } else {
            return ( R::pnorm( 0.0, xvals.val_max, bandwidth, 1, 0 ) );
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
// [[Rcpp::export]]
double ProbabilityLessSmoothed(
    const Rcpp::NumericVector&  X,
    const Rcpp::NumericVector&  Y,
    int   nsteps = 100,
    double bandwidth = NA_REAL
){
    if ( X.size() == 0 ) {
        throw Rcpp::exception( "X is empty" );
    }
    if ( Y.size() == 0 ) {
        throw Rcpp::exception( "Y is empty" );
    }
    ImportedValues xvals( X );
    ImportedValues yvals( Y );
    BinnedValues diffBins = BinnedValues::difference( xvals, yvals, nsteps );

    return ( diffBins.probabilityNonPositive( bandwidth ) );
}

typedef std::map<std::string, std::size_t> name_index_map_t;

name_index_map_t nameToIndexMap( const Rcpp::CharacterVector& strings )
{
    name_index_map_t res;
    for ( std::size_t i = 0; i < (std::size_t)strings.size(); i++ ) {
        res.insert( std::make_pair( Rcpp::as<std::string>( strings[i] ), i ) );
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
// [[Rcpp::export]]
Rcpp::List DifferenceStatistics(
    Rcpp::NumericMatrix&  X,
    Rcpp::NumericMatrix&  Y,
    const Rcpp::NumericVector&  Deltas,
    int   nsteps = 100,
    double maxBandwidth = NA_REAL,
    bool  matchIterations = true
){
    if ( X.ncol() == 0 || X.nrow() == 0 ) {
        throw Rcpp::exception( "X matrix is empty" );
    }
    if ( Y.ncol() == 0 || Y.nrow() == 0 ) {
        throw Rcpp::exception( "Y matrix is empty" );
    }
    if ( !R_IsNA( maxBandwidth ) && (maxBandwidth < 0.0) ) {
        throw Rcpp::exception( "maxBandwidth cannot be negative" );
    }

    Rcpp::Dimension deltaDims = Deltas.attr("dim");
    if ( deltaDims.size() != 3 ) {
        throw Rcpp::exception( "Deltas is not 3D array" );
    }
    if ( deltaDims[2] != X.nrow() ) {
        throw Rcpp::exception( "Number of rows in X and 3rd dimension in Deltas does not match" );
    }
    if ( deltaDims[2] != Y.nrow() ) {
        throw Rcpp::exception( "Number of rows in Y and 3rd dimension in Deltas does not match" );
    }
    Rcpp::List deltaDimNames = Deltas.attr("dimnames");
    LOG_DEBUG2( "Loading delta 1st dimension names" );
    name_index_map_t deltaXNames = nameToIndexMap( deltaDimNames[0] );
    LOG_DEBUG2( "Loading delta 2nd dimension names" );
    name_index_map_t deltaYNames = nameToIndexMap( deltaDimNames[1] );
    std::size_t deltaStride = deltaDims[0] * deltaDims[1]; // stride of vectors in the 3rd dimension
    std::size_t nsamples = deltaDims[2];

    LOG_DEBUG2( "Loading X column names" );
    name_vec_t xColNames = Rcpp::as<name_vec_t>( Rcpp::as<Rcpp::List>( X.attr("dimnames") )[1] );
    LOG_DEBUG2( "Loading Y column names" );
    name_vec_t yColNames = Rcpp::as<name_vec_t>( Rcpp::as<Rcpp::List>( Y.attr("dimnames") )[1] );

    std::vector<double> xvals_merged;
    std::vector<double> yvals_adj_merged;
    xvals_merged.reserve( X.ncol() * Y.ncol() * nsamples );
    yvals_adj_merged.reserve( X.ncol() * Y.ncol() * nsamples );

    for ( int i = 0; i < X.ncol(); i++ ) {
        ImportedValues xvals( X.column(i) );

        name_index_map_t::const_iterator iColIt = deltaXNames.find( xColNames[i] );
        if ( iColIt == deltaXNames.end() ) {
            THROW_EXCEPTION( Rcpp::exception,
                             "Column '" << xColNames[i] << "' not found in Deltas 1st dimension" );
        }

        for ( int j = Y.ncol()-1; j >= 0; j-- ) {
            LOG_DEBUG2( "Processing distribution of " << xColNames[i] << " vs " << yColNames[j] );
            name_index_map_t::const_iterator jColIt = deltaYNames.find( yColNames[j] );
            if ( jColIt == deltaYNames.end() ) {
                THROW_EXCEPTION( Rcpp::exception,
                                 "Column '" << yColNames[j] << "' not found in Deltas 2nd dimension" );
            }
            xvals_merged.insert( xvals_merged.end(), xvals.values.begin(), xvals.values.end() );

            Rcpp::MatrixColumn<REALSXP> yCol = Y.column(j);
            std::size_t deltaOffset = iColIt->second + jColIt->second * deltaDims[0];
            for ( std::size_t k = 0; k < nsamples; k++ ) {
                yvals_adj_merged.push_back( yCol[k] + Deltas[ deltaOffset + k * deltaStride ] );
            }
        }
    }
    LOG_DEBUG2( "Calculating probability that difference is non-positive" );
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
        if ( !R_IsNA( maxBandwidth ) && (bw > maxBandwidth) ) {
            bw = maxBandwidth;
        }
        prob = diffBins.probabilityNonPositive( bw );
        diff_mean = diffBins.average();
        diff_var = diffBins.variance();
    } else {
        // bandwidth is calculated from the global difference of X values vs Y values
        BinnedValues globalDiffBins = BinnedValues::difference( ImportedValues( xvals_merged ),
                                                                ImportedValues( yvals_adj_merged ), nsteps );
        bw = maxBandwidth == 0.0 ? 0.0 : globalDiffBins.bw_nrd();
        if ( !R_IsNA( maxBandwidth ) && (bw > maxBandwidth) ) {
            bw = maxBandwidth;
        }
        // calculate probability for every X column and all Y columns
        // NOTE it's not correct to calculate the probability using globalDiffBins,
        // because every X-Y column pair has its own adjustment
        std::size_t chunk_size = Y.ncol() * nsamples;
        std::vector<double> probs;
        double max_prob = 0;
        double sum_mean = 0;
        double sum_var = 0;
        for ( int i = 0; i < X.ncol(); i++ ) {
            BinnedValues diffBins = BinnedValues::difference( ImportedValues( std::vector<double>( xvals_merged.begin() + i * chunk_size,
                                                                              xvals_merged.begin() + (i+1) * chunk_size ) ),
                                                 ImportedValues( std::vector<double>( yvals_adj_merged.begin() + i * chunk_size,
                                                                              yvals_adj_merged.begin() + (i+1) * chunk_size ) ), nsteps );
            probs.push_back( diffBins.probabilityNonPositive( bw ) );
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
    LOG_DEBUG( "P=" << prob );
    return ( Rcpp::List::create(
        Rcpp::Named( "mean", diff_mean ),
        Rcpp::Named( "sd", sqrt( diff_var ) ),
        Rcpp::Named( "prob_less", prob ),
        Rcpp::Named( "bw", bw )
        ) );
}

//??? Calculates average probabilities that given contrasts would be
//??? be less or equal than zero.
//??? Gaussian kernel smoothing is used to represent the contrast distribution.
//??? The probability is averaged across all possible combinations of
//??? the experiments of the relevant conditions.
//???
//??? @experiments2conditions_experiment index of experiment in many-to-many map
//??? @experiments2conditions_condition index of condition in many-to-many map
//??? @contrastXcondition contrast matrix, rows are contrasts, columns are conditions
//??? @X matrix of samples of X random variable, columns are different experiments, rows are MCMC iterations
//??? @nsteps the number of segments to divide the range of X values into
//??? @maxBandwidth the maximum gaussian smoothing kernel bandwidth, 0 disables smoothing
//??? @return E(X), D(X), P(X>=0), P(X<=0)
// [[Rcpp::export]]
Rcpp::List ContrastStatistics(
    const Rcpp::NumericMatrix&  X,
    const Rcpp::IntegerVector&  experiments2conditions_experiment,
    const Rcpp::IntegerVector&  experiments2conditions_condition,
    const Rcpp::NumericMatrix&  contrastXcondition,
    int   nsteps = 100,
    double maxBandwidth = NA_REAL,
    const Rcpp::NumericVector&  quant_probs = Rcpp::NumericVector::create(0.025, 0.25, 0.50, 0.75, 0.975)
){
    std::size_t ncontrasts = contrastXcondition.nrow();
    std::size_t nconditions = contrastXcondition.ncol();
    LOG_DEBUG2("ncontrasts=" << ncontrasts << " nconditions=" << nconditions);
    if ( contrastXcondition.ncol() == 0 || contrastXcondition.nrow() == 0 ) {
        THROW_EXCEPTION( std::invalid_argument, "contrastXcondition matrix is empty" );
    }
    if ( X.ncol() == 0 || X.nrow() == 0 ) {
        THROW_EXCEPTION( std::invalid_argument, "X matrix is empty" );
    }
    if ( experiments2conditions_experiment.size() != experiments2conditions_condition.size() ) {
        THROW_EXCEPTION(std::invalid_argument, "experiments2conditions: experiment and condition vectors have different number of elements" );
    }
    if ( !R_IsNA( maxBandwidth ) && (maxBandwidth < 0.0) ) {
        THROW_EXCEPTION(std::invalid_argument, "maxBandwidth cannot be negative" );
    }

    std::size_t nsamples = X.nrow();
    std::size_t nexperiments = X.ncol();
    LOG_DEBUG2("nmcmc_samples=" << nsamples << " nexperiments=" << nexperiments);

    experiment_map_t condition2experiments;

    for ( int i = 0; i < experiments2conditions_condition.size(); i++ ) {
        std::size_t cond_ix = experiments2conditions_condition[i];
        std::size_t exper_ix = experiments2conditions_experiment[i];
        if ( cond_ix == 0 ) {
            THROW_EXCEPTION(std::invalid_argument, "Condition index should be positive and <=" << nconditions <<", " << cond_ix << " found");
        } else {
            --cond_ix;
        }
        if ( exper_ix == 0 || exper_ix > nexperiments ) {
          THROW_EXCEPTION(std::invalid_argument, "Experiment index should be positive and <=" << nexperiments << ", " << exper_ix << " found");
        } else {
          --exper_ix;
        }
        experiment_map_t::iterator condIt = condition2experiments.find( cond_ix );
        if ( condIt == condition2experiments.end() ) {
            condition2experiments.insert( condIt, std::make_pair( cond_ix, experiment_set_t(1, exper_ix) ) );
        } else {
            condIt->second.push_back(exper_ix);
        }
    }

    std::vector<std::size_t> index_contrast(ncontrasts, 0);
    std::vector<double> prob_nonpos(ncontrasts, NAN);
    std::vector<double> prob_nonneg(ncontrasts, NAN);
    std::vector<double> contrast_mean(ncontrasts, NAN);
    std::vector<double> contrast_sd(ncontrasts, NAN);
    std::vector<double> bw(ncontrasts, NAN);
    std::vector<double> bin_width(ncontrasts, NAN);
    std::vector<std::size_t> ncombn(ncontrasts, 0);
    std::vector<std::size_t> nsamp(ncontrasts, 0);
    std::vector<std::vector<double> > quantiles(quant_probs.size(), std::vector<double>(ncontrasts, NAN));

    // generate constrast samples for all relevant experiments combinations
    std::vector<double> all_contrast_samples;
    std::vector<double> experiment_contrast_samples(nsamples, 0.0);

    LOG_DEBUG2( "Calculating the probabilities that contrast is non-positive and non-negative" );
    for ( std::size_t contr_ix = 0; contr_ix < ncontrasts; ++contr_ix ) {
        // calculate how many conditions are associated with current contrast
        std::vector<std::size_t> cond_ixs;
        std::size_t cur_ncombn = 1;
        for ( std::size_t cond_ix = 0; cond_ix < nconditions; ++cond_ix ) {
            if (contrastXcondition.at(contr_ix, cond_ix) != 0.0) {
                LOG_DEBUG2("contrastXcondition["<<contr_ix<<", "<<cond_ix<<"]="<<contrastXcondition.at(contr_ix, cond_ix));
                cond_ixs.push_back( cond_ix );
                experiment_map_t::const_iterator condIt = condition2experiments.find(cond_ix);
                if ( condIt != condition2experiments.end() ) {
                    cur_ncombn *= condIt->second.size();
                } else {
                    THROW_EXCEPTION(std::invalid_argument, "No experiments for condition #" << (cond_ix+1) << " found");
                }
            }
        }
        if (cur_ncombn == 0) { // should not happen as checked before
            THROW_EXCEPTION(Rcpp::exception, "No experiment combinations for contrast #" << contr_ix);
        }
        std::vector<std::size_t> next_exper_i(cond_ixs.size(), 0); // counter of experiments for each condition
        all_contrast_samples.resize(0);
        all_contrast_samples.reserve(nsamples * cur_ncombn);
        for (std::size_t comb_i = 0; comb_i < cur_ncombn; ++comb_i) {
            // generate contrast samples for all experiment combinations
            // relevant to a given contrast
            std::fill(experiment_contrast_samples.begin(), experiment_contrast_samples.end(), 0.0);
            for ( std::size_t cond_i = 0; cond_i < cond_ixs.size(); ++cond_i ) {
                std::size_t cond_ix = cond_ixs[cond_i];
                std::size_t exper_ix = condition2experiments.find(cond_ix)->second[next_exper_i[cond_i]];
                double cond_w = contrastXcondition.at(contr_ix, cond_ix);
                for ( std::size_t i  = 0; i < nsamples; ++i ) {
                    experiment_contrast_samples[i] += cond_w * X.at(i, exper_ix);
                }
            }
            // push the contrasts for current experiments into global vector
            all_contrast_samples.insert(all_contrast_samples.end(),
                experiment_contrast_samples.begin(), experiment_contrast_samples.end());
            // increment experiments enumeration
            for ( std::size_t cond_i = 0; cond_i < cond_ixs.size(); ++cond_i ) {
                if (next_exper_i[cond_i] + 1 < condition2experiments.find(cond_ixs[cond_i])->second.size()) {
                    next_exper_i[cond_i]++;
                    break;
                } else {
                    next_exper_i[cond_i] = 0;
                }
            }
        }
        // calculate the quantiles
        bacc::accumulator_set<double, bacc::stats<bacc::tag::extended_p_square> >
            quant_acc( bacc::tag::extended_p_square::probabilities = Rcpp::as<std::vector<double> >(quant_probs) );
        for ( std::size_t i = 0; i < all_contrast_samples.size(); ++i ) {
            quant_acc( all_contrast_samples[i] );
        }

        // calculate the contrast
        ImportedValues contrast_vals( all_contrast_samples );
        BinnedValues contrast_bins( contrast_vals, contrast_vals.defaultBinWidth( nsteps ) );

        double cur_bw = maxBandwidth == 0.0 ? 0.0 : contrast_bins.bw_nrd();
        if ( !R_IsNA( maxBandwidth ) && (cur_bw > maxBandwidth) ) {
            cur_bw = maxBandwidth;
        }
        index_contrast[contr_ix] = contr_ix + 1;
        prob_nonpos[contr_ix] = contrast_bins.probabilityNonPositive( cur_bw );
        prob_nonneg[contr_ix] = contrast_bins.probabilityNonNegative( cur_bw );
        contrast_mean[contr_ix] = contrast_bins.average();
        contrast_sd[contr_ix] = sqrt(contrast_bins.variance());
        bw[contr_ix] = cur_bw;
        bin_width[contr_ix] = contrast_bins.step;
        ncombn[contr_ix] = cur_ncombn;
        nsamp[contr_ix] = all_contrast_samples.size();
        for ( int i = 0; i < quant_probs.size(); ++i ) {
            quantiles[i][contr_ix] = bacc::extended_p_square(quant_acc)[i];
        }
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
        LOG_DEBUG2(contr_ix << "-th contrast done");
    }
    Rcpp::List res = Rcpp::List::create(
        Rcpp::Named( "index_contrast", index_contrast ),
        Rcpp::Named( "mean", contrast_mean ),
        Rcpp::Named( "sd", contrast_sd ),
        Rcpp::Named( "prob_nonneg", prob_nonneg ),
        Rcpp::Named( "prob_nonpos", prob_nonpos ),
        Rcpp::Named( "bin_width", bin_width ),
        Rcpp::Named( "bw", bw ),
        Rcpp::Named( "ncombn", ncombn ),
        Rcpp::Named( "nsamples", nsamp )
        );
    for ( int i = 0; i < quant_probs.size(); ++i ) {
        std::stringstream ss;
        ss << quant_probs[i]*100 << "%";
        res[ss.str()] = quantiles[i];
    }
    return res;
}
