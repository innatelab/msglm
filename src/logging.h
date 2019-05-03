#pragma once

#include <sstream>
#include <Rcpp.h>

#ifndef NDEBUG
//#define DYNLOAD_DEBUG
#endif

#define THROW_EXCEPTION( excp_class, msg ) { \
    std::ostringstream __excp__msg__; \
    __excp__msg__ << msg; \
    throw excp_class( __excp__msg__.str().c_str() ); \
}

#define LOG_LEVEL 0

#define LOG_RCOUT( msg ) Rcpp::Rcout << msg << "\n"

#if LOG_LEVEL >= 0
#define LOG_DEBUG0( msg ) LOG_RCOUT( msg )
#else
#define LOG_DEBUG0( msg )
#endif

#if LOG_LEVEL >= 1
#define LOG_DEBUG1( msg ) LOG_RCOUT( msg )
#else
#define LOG_DEBUG1( msg )
#endif
#define LOG_DEBUG(msg) LOG_DEBUG1(msg)

#if LOG_LEVEL >= 2
#define LOG_DEBUG2( msg ) LOG_RCOUT( msg )
#else
#define LOG_DEBUG2( msg )
#endif
