#pragma once

#include <cmath>
#include <Rcpp.h>

#include "imported_values.h"

// counts of observations distributed to bins
typedef std::vector<std::size_t> bins_t;

// vector of dimension names
typedef std::vector<std::string> name_vec_t;

// distributes the values of the input vector
// to the bins according to the specified bin width and rounding mode
struct BinnedValues {
    double          val_min;
    double          val_max;

    double          step;
    size_t          bins_sum;
    bins_t          bins;

    BinnedValues(const ImportedValues& values, double binWidth, bool floor = true, double max_support_expand = 10.0);

    // uninitialized
    BinnedValues()
    : val_min( std::numeric_limits<value_t>::quiet_NaN() )
    , val_max( std::numeric_limits<value_t>::quiet_NaN() )
    , step( std::numeric_limits<value_t>::quiet_NaN() )
    , bins_sum( 0 )
    {
    }

    // degenerated
    BinnedValues(double val)
    : val_min( val ), val_max( val ), bins_sum( 1 ), bins( 1, 1 )
    {
    }

    BinnedValues(const bins_t& bins, double binWidth, double val_min, double val_max);

    std::size_t size() const { return bins.size(); }

    // Probability that random variable
    // would be less or equal than zero using the
    // Gaussian kernel smoothing.
    // @bins binned samples of random variable
    // @start value corresponding to the first bin
    // @step value step between the bins
    // @bandwidth the Gaussian smoothing kernel bandwidth, defaults to the segment size, 0 disables smoothing
    // @return P(X<=0) if negative, P(X>=0) if !negative
    double compareWithZero(double bandwidth = NA_REAL, bool negative = true) const;

    double probabilityNonPositive(double bandwidth = NA_REAL) const {
        return compareWithZero(bandwidth, true);
    }
    double probabilityNonNegative(double bandwidth = NA_REAL) const {
        return compareWithZero(bandwidth, false);
    }

    double norm_average() const;
    double average() const {
        return bins_sum > 0
                ? val_min + step * norm_average()
                : 0.5 * (val_max + val_min);
    }

    // returns the variance normalized by step
    // to get the variance, multiply it by step^2
    double norm_variance() const;
    double variance() const { return step*step * norm_variance(); }
    double sd() const { return step * sqrt(norm_variance()); }

    // rule-of-thumb method for bandwidth selection
    // returns bw normalized by step, to get the bw, multiply it by step
    // see R help: bw.nrd()
    static double norm_bw_nrd(double bins_sum, double sd, double quartile1, double quartile3);

    // rule-of-thumb method for bandwidth selection
    // returns bw normalized by step, to get the bw, multiply it by step
    // see R help: bw.nrd()
    double norm_bw_nrd() const;

    // rule-of-thumb method for bandwidth selection
    // see R help: bw.nrd()
    double bw_nrd() const { return step * norm_bw_nrd(); }

    // make bins for the distribution of X-Y
    static BinnedValues difference(const ImportedValues& xvals,
                                   const ImportedValues& yvals,
                                   std::size_t nsteps
    );
};
