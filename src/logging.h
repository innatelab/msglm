#pragma once

#include <sstream>
#include <cpp11/R.hpp>
#include <cpp11/function.hpp>

#ifndef NDEBUG
//#define DYNLOAD_DEBUG
#endif

#define THROW_EXCEPTION(excp_class, msg, ...) { \
    char __excp_buff__[1024];                   \
    int __excp_msg_len__ = std::snprintf(__excp_buff__, 1024, msg, ##__VA_ARGS__); \
    if (__excp_msg_len__ >= 0 && __excp_msg_len__ < 1024) { \
        throw excp_class(__excp_buff__);        \
    } else {                                    \
        throw excp_class(msg);                  \
    }                                           \
}

#define LOG_LEVEL 0

#define LOG_RCOUT(msg, ...) cpp11::message(msg, ##__VA_ARGS__)

#if LOG_LEVEL >= 0
#define LOG_DEBUG0(msg, ...) LOG_RCOUT(msg, __VA_ARGS__)
#else
#define LOG_DEBUG0(msg, ...)
#endif

#if LOG_LEVEL >= 1
#define LOG_DEBUG1(msg, ...) LOG_RCOUT(msg, __VA_ARGS__)
#else
#define LOG_DEBUG1(msg, ...)
#endif
#define LOG_DEBUG(msg, ...) LOG_DEBUG1(msg, __VA_ARGS__)

#if LOG_LEVEL >= 2
#define LOG_DEBUG2(msg, ...) LOG_RCOUT(msg, __VA_ARGS__)
#else
#define LOG_DEBUG2(msg, ...)
#endif
