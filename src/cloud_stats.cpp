// [[Rcpp::depends(BH)]]

#include "binned_values.h"

#include "logging.h"

#include <algorithm>
#include <cmath>
#include <boost/math/distributions/normal.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/extended_p_square.hpp>

#include <cpp11.hpp>

using boost::math::normal; // typedef provides default type is double.
using namespace cpp11;

#if 0 // unsupported default params
writable::data_frame CloudsOverlap(
        doubles_matrix<by_column>  X,
        doubles_matrix<by_column>  Y,
        double maxBandwidth = na<double>(),
        int   nsteps = 100,
        doubles quant_probs = writable::doubles{0.025, 0.25, 0.50, 0.75, 0.975}
)
#endif

[[cpp11::register]]
writable::data_frame CloudsOverlap(
    doubles_matrix<by_column>  X,
    doubles_matrix<by_column>  Y,
    double maxBandwidth,
    int   nsteps,
    doubles quant_probs
){
    //LOG_DEBUG0("CloudsOverlap");
    if (X.ncol() == 0 || X.nrow() == 0) {
        THROW_EXCEPTION(std::invalid_argument, "X matrix is empty");
    }
    if (Y.ncol() == 0 || Y.nrow() == 0) {
        THROW_EXCEPTION(std::invalid_argument, "Y matrix is empty");
    }
    if (X.ncol() != Y.ncol() || X.nrow() != Y.nrow()) {
        THROW_EXCEPTION(std::invalid_argument, "Y matrix dimensions (%d, %d) are different from X (%d, %d)",
                        Y.nrow(), Y.ncol(), X.nrow(), X.ncol());
    }
    if (!R_IsNA(maxBandwidth) && (maxBandwidth < 0.0)) {
        THROW_EXCEPTION(std::invalid_argument, "maxBandwidth cannot be negative");
    }

    const std::size_t ndims = X.ncol();
    const std::size_t nsamples = X.nrow();

    // calculate bandwidth for each X dimension (independently)
    std::vector<double> X_bw(ndims);
    std::vector<double> X_bw_inv(ndims);
    for (std::size_t k = 0; k < ndims; ++k) {
        const ImportedValues X_col(X[k]);
        const BinnedValues X_bins(X_col, X_col.defaultBinWidth(nsteps));
        double bw = maxBandwidth == 0.0 ? 0.0 : X_bins.bw_nrd();
        if (!is_na(maxBandwidth) && (bw > maxBandwidth)) {
            bw = maxBandwidth;
        }
        X_bw[k] = bw;
        X_bw_inv[k] = 1.0/bw;
        //LOG_DEBUG0("X_bw[" << k << "]=" << bw);
    }

    // check every X-Y pair
    double prob_sum = 0.0;
    double offset_sum = 0.0;
    double scale_sum = 0.0;
    std::vector<double> delta(ndims);
    for (std::size_t i = 0; i < nsamples; ++i) {
        //LOG_DEBUG0("i=%d", i);
        // the direction to the X set
        double scale = 0.0;
        for (std::size_t k = 0; k < ndims; ++k) {
            //LOG_DEBUG0("X[%d, %d]=%g Y[%d, %d]=%g", i, k, X(i, k), i, k, Y(i, k));
            delta[k] = (Y(i,k) - X(i,k))*X_bw_inv[k];
            scale += delta[k]*delta[k];
            delta[k] *= X_bw_inv[k];
        }
        scale = sqrt(scale);
        if (scale > 0.0) {
            //LOG_DEBUG0("scale=%g", scale);
            // loop over all points of X and calculate the probability that
            // random point generated around Xj would be further than Yi
            double Y_further_prob_sum = 0.0;
            double i_offset_sum = 0.0;
            boost::math::normal distr(0.0, scale);
            for (std::size_t j = 0; j < nsamples; ++j) {
                double d = 0.0;
                for (std::size_t k = 0; k < ndims; ++k) {
                    d += delta[k]*(Y(i,k)-X(j,k));
                }
                i_offset_sum += d;
                Y_further_prob_sum += cdf(complement(distr, d));//d > 0 ? 0.0 : 1.0;
            }
            scale_sum += scale;
            offset_sum += i_offset_sum/nsamples;
            prob_sum += Y_further_prob_sum/nsamples;
        } else {
            prob_sum += 0.5/(ndims+2);
        }
    }
    //LOG_DEBUG0("prob_avg=%g", prob_sum/nsamples);
    //LOG_DEBUG0("offset_avg=%g", offset_sum/nsamples);
    //LOG_DEBUG0("scale_avg=%g", scale_sum/nsamples);
    return ( writable::data_frame{
              "offset"_nm = writable::doubles{offset_sum/nsamples},
              "scale"_nm = writable::doubles{scale_sum/nsamples},
              "prob_further"_nm = writable::doubles{prob_sum/nsamples} } );
}
