#include "binned_values.h"

#include <cpp11.hpp>
#include <algorithm>
#include <numeric>
#include <boost/accumulators/statistics/stats.hpp>
#include "logging.h"

#include "Rmath.h"

using namespace cpp11;

BinnedValues::BinnedValues(const bins_t& bins, double binWidth, double val_min, double val_max)
    : val_min( val_min ), val_max( val_max )
    , step( binWidth )
    , bins_sum( std::accumulate( bins.begin(), bins.end(), 0 ) )
    , bins( bins )
{
}
BinnedValues::BinnedValues(const ImportedValues& values, double binWidth,
                           bool floor, double max_support_expand)
    // define the bin bounds, avoid overexpanding the values support (in case of very distant outliers)
    : val_min( std::max(values.val_min, 0.5*(values.support_min+values.support_max) - 0.5*max_support_expand*(values.support_max-values.support_min)) )
    , val_max( std::min(values.val_max, 0.5*(values.support_min+values.support_max) + 0.5*max_support_expand*(values.support_max-values.support_min)) )
    , step( binWidth )
    , bins_sum( values.values.size() )
{
    if ( binWidth < 0 ) THROW_EXCEPTION(std::invalid_argument, "step=%g is negative", binWidth);
    if ( val_min == val_max ) {
        // degenerated case, everything concentrated at single value
        bins.resize( 1, 0 );
        bins[0] = bins_sum;
        return;
    }
    else {
        if ( binWidth == 0 ) THROW_EXCEPTION(std::invalid_argument, "step is zero");
        bins.resize( floor ? ( val_max - val_min ) / binWidth + 1
                     : std::ceil( ( val_max - val_min ) / binWidth ) + 1, 0 );
    }
    for ( std::size_t i = 0; i < values.values.size(); ++i ) {
        const double val = values.values[i];
        // quantize val and project it into [0, bins.size()-1] range
        int ix = std::min(std::max((int)(floor
                                         ? (val - val_min)/step
                                         : std::ceil((val_max - val)/step)), 0), (int)bins.size()-1);
#if 0
        if ( ix < 0 || ix >= bins.size() ) {
            THROW_EXCEPTION(std::out_of_range,
                            "Bin #%d for value %g is out of bounds [%g, %g] (nbins=%d)",
                            i, val, val_min, val_max, bins.size());
        }
#endif
        bins[ix]++;
    }
#if 0
    std::size_t bins_sum_ = std::accumulate( bins.begin(), bins.end(), 0 );
    LOG_DEBUG2( "Sum is " << bins_sum_ );
    if ( bins_sum_ != bins_sum ) {
        THROW_EXCEPTION(std::length_error, "%d element(s) in bins, %d expected",
                        bins_sum_, bins_sum);
    }
#endif
}

// Probability that random variable
// would be less or equal than the given value using the
// Gaussian kernel smoothing.
// @y value to compare with
// @bandwidth the Gaussian smoothing kernel bandwidth, defaults to the segment size, 0 disables smoothing
// @return P(X<=y) if negative, P(X>=y) if !negative
double BinnedValues::probabilityCompareWith(double y, double bandwidth, bool negative) const
{
    LOG_DEBUG2("minmax=[%g, %g] step=%g size=%d", val_min, val_max, step, size())

    if ( val_max == val_min ) {
        // distribution is degenerated
        LOG_DEBUG1("Degenerated distribution with minmax=[%g, %g]", val_min, val_max);
        return ( is_na( bandwidth )
                ? ( val_max <= y ? 1.0 : 0.0 )
                : Rf_pnorm5( y, val_max, bandwidth, negative, 0 ) );
    }
    else if ( (negative && ( ( val_max - y < -5.0*step*size() ) || ( val_min - y > 30.0*step*size() ) ))
           || (!negative && ( ( val_max - y < -30.0*step*size() ) || ( val_min - y > 5.0*step*size() ) ))
    ){
        // distribution is almost degenerated w.r.t. distance to zero
        // (30 times the value range)
        LOG_DEBUG1("Almost degenerated distribution");
        if ( is_na( bandwidth ) ) {
            bandwidth = sqrt( norm_variance() );
        } else {
            bandwidth /= step;
        }
        LOG_DEBUG1("Normalized bandwidth=%g", bandwidth);
        return ( Rf_pnorm5( 0.0, (val_min - y + norm_average()) / step, bandwidth, negative, 0 ) );
    }
    if ( bins_sum == 0 ) return ( 0.5 );

    // if bandwidth is not specified, use the rule-of-thumb
    if ( is_na( bandwidth ) ) bandwidth = norm_bw_nrd();
    else {
        if ( bandwidth < 0 ) throw std::invalid_argument("Negative bandwidth not allowed");
        bandwidth /= step;
    }
    LOG_DEBUG2("Normalized bandwidth=%g", bandwidth);

    double res = 0.0;
    const double offset = (val_min - y) / step + 0.5;
    if ( bandwidth > 0.0 ) {
        // integrate the Gaussian kernel probability across all bins
        for ( std::size_t i = 0; i < bins.size(); ++i ) {
            const size_t bin_i = bins[i];
            if (bin_i > 0) res += bin_i * Rf_pnorm5( 0.0, offset + i, bandwidth, negative, 0 );
        }
    } else {
        // no kernel, count bins corresponding to non-positive differences
        for ( std::size_t i = std::max(0, (int)(-offset)-1); i < bins.size(); ++i ) {
            if ( offset + i <= 0.0 ) res += bins[i];
        }
        if ( !negative ) res = bins_sum - res;
    }
    LOG_DEBUG2("P(X %c= %g) = %g", negative ? '<' : '>', y, res / bins_sum);
    return res / bins_sum;
}

double BinnedValues::norm_average() const {
    if ( bins_sum > 0 ) {
        double res = 0.0;
        for ( std::size_t i = 0; i < bins.size(); ++i ) {
            res += bins[i] * i;
        }
        return res / bins_sum + 0.5;
    } else {
        return 0.0;
    }
}

double BinnedValues::norm_variance() const {
    if ( bins_sum > 0 ) {
        double res = 0.0;
        double avg = norm_average() - 0.5;
        for ( std::size_t i = 0; i < bins.size(); ++i ) {
            res += bins[i] * (i-avg) * (i-avg);
        }
        return res / bins_sum;
    } else {
        return 0.0;
    }
}

// rule-of-thumb method for bandwidth selection
// returns bw normalized by step, to get the bw, multiply it by step
// see R help: bw.nrd()
double BinnedValues::norm_bw_nrd(double bins_sum, double sd, double quartile1, double quartile3)
{
    const double h = (quartile3 > quartile1 ? quartile3 - quartile1 : 1) / 1.34;
    const double a = sd < h ? sd : h;
    const double res = a / pow(0.75 * bins_sum, 0.2);
    LOG_DEBUG2("interquartile=(%g, %g) h=%g sd=%g bw=%g", quartile1, quartile3, h, sd, res);
    return res;
}

// rule-of-thumb method for bandwidth selection
// returns bw normalized by step, to get the bw, multiply it by step
// see R help: bw.nrd()
double BinnedValues::norm_bw_nrd() const
{
    size_t i_quartile_1 = 0;
    size_t i_quartile_3 = 0;
    size_t n_elems = 0;
    for ( std::size_t i = 0; i < bins.size(); ++i ) {
        const size_t n_new_elems = n_elems + bins[i];
        if ( (4*n_elems < bins_sum) && (4*n_new_elems >= bins_sum) ) {
            i_quartile_1 = i;
        }
        if ( (4*n_elems < 3*bins_sum) && (4*n_new_elems >= 3*bins_sum) ) {
            i_quartile_3 = i;
            break;
        }
        n_elems = n_new_elems;
    }
    return norm_bw_nrd(bins_sum, sqrt(norm_variance()), i_quartile_1, i_quartile_3);
}

// make bins for the distribution of X-Y
BinnedValues BinnedValues::difference( const ImportedValues& xvals,
                                       const ImportedValues& yvals,
                                       std::size_t nsteps
){
    double step = ( ( xvals.val_max - xvals.val_min )
                    + ( yvals.val_max - yvals.val_min ) ) / nsteps;

    BinnedValues xbins( xvals, step, true );
    BinnedValues ybins( yvals, step, false );

    // use the bins of X and Y to calculate the distribution of X-Y
    bins_t diff( xbins.size() + ybins.size(), 0 );
    for ( std::size_t i = 0; i < xbins.size(); ++i ) {
        std::size_t x = xbins.bins[i];
        if ( x > 0 ) for ( std::size_t j = 0; j < ybins.size(); ++j ) {
            std::size_t y = ybins.bins[j];
            if ( y > 0 ) {
#if 0
                if ( i + j >= diff.size() ) {
                    THROW_EXCEPTION(std::out_of_range,
                                    "Diff.Bin index %d+%d=%d is out of bounds (%d)",
                                    i, j, i+j, diff.size());
                }
#endif
                diff[ i + j ] += x * y;
            }
        }
    }
    return ( BinnedValues( diff, step, xbins.val_min - ybins.val_max,
                           xbins.val_max - ybins.val_min ) );
}
