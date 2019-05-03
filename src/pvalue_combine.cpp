#include <algorithm>
#include <map>
#include <Rcpp.h>
// [[Rcpp::depends(BH)]]
#include <boost/math/tools/roots.hpp>
#include <Rcpp/Rmath.h>

#ifndef NDEBUG
//#define DYNLOAD_DEBUG
#endif

typedef std::map<std::size_t, double> inner_param_map_t;
typedef std::map<double, inner_param_map_t> param_map_t;

static param_map_t GammaMethodParams;

class ParamEquation {
    double      stt;
    std::size_t nTests;

public:
    ParamEquation(
        double      stt,
        std::size_t nTests
    ) : stt( stt ), nTests( nTests )
    {}

    double operator()( double x ) const
    {
        return ( R::pgamma( nTests * R::qgamma( stt, x, 1.0, 0, 0 ),
                            nTests * x, 1.0, 0, 0 ) - stt );
    }
};

double GammaMethodParam(
    double      stt,
    std::size_t nTests
){
    // find param in cache
    param_map_t::iterator pmit = GammaMethodParams.find( stt );
    if ( pmit == GammaMethodParams.end() ) {
        pmit = GammaMethodParams.insert( pmit, std::make_pair( stt, inner_param_map_t() ) );
    }
    inner_param_map_t::iterator ipmit = pmit->second.find( nTests );
    if ( ipmit != pmit->second.end() ) {
        return ( ipmit->second );
    }
    boost::uintmax_t maxIter = 10000;
    std::pair<double, double> rootBracket = boost::math::tools::toms748_solve(
        ParamEquation( stt, nTests ), 1E-15, 10.0, boost::math::tools::eps_tolerance<double>(20), maxIter );
    double res = 0.5 * ( rootBracket.first + rootBracket.second );
    pmit->second.insert( ipmit, std::make_pair( nTests, res ) );
    return ( res );
}
    
//??? Combine p-values according to Zaykin (2007) Pharm.Stat. (10.1002/pst.304).
//???
//??? @param pvalues P-values to combine
//??? @param stt Soft Truncation Threshold
//??? @return combined p-value
// [[Rcpp::export]]
double GammaMethodPValuesCombine( std::vector<double> pvalues, double stt )
{
    if ( pvalues.size() == 1 ) return ( pvalues[0] );
    double param = GammaMethodParam( stt, pvalues.size() );
    double Y = 0;
    for ( std::size_t i = 0; i < pvalues.size(); ++i ) {
        Y += R::qgamma( pvalues[i], param, 1.0, 0, 0 );
    }
    return ( R::pgamma( Y, param * pvalues.size(), 1.0, 0, 0 ) );
}

#if 0
//??? Calculate parameter for GammaMethod of p-values combining.
//???
//??? @param nTests number of tests
//??? @param stt Soft Truncation Threshold
//??? @return parameter
/// [[Rcpp::export]]
double GammaMethodPValuesCombineParam( std::size_t nTests, double stt )
{
    return ( GammaMethodParam( stt, nTests ) );
}
#endif