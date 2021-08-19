#pragma once

#include <sstream>
#include <cpp11/R.hpp>
#include <cpp11/function.hpp>

#ifndef NDEBUG
//#define DYNLOAD_DEBUG
#endif

#define THROW_EXCEPTION( excp_class, msg ) { \
    std::ostringstream __excp__msg__; \
    __excp__msg__ << msg; \
    throw excp_class( __excp__msg__.str() ); \
}

#define LOG_LEVEL 0

#define LOG_RCOUT( msg ) { \
  std::ostringstream __msg__;       \
  __msg__ << msg;                   \
  cpp11::message(__msg__.str()); \
}

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
