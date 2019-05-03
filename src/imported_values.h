#pragma once

#include <vector>
#include <Rcpp.h>

// values
typedef double value_t;
typedef std::vector<value_t> values_t;

// stores the vector and its min/max values
// serves an input to BinnedValues
struct ImportedValues {
    const values_t  values;
    double          val_min;
    double          val_max;
    double          support_min;
    double          support_max;

    void init_bounds(double support_prob);

    ImportedValues(const Rcpp::NumericVector& values, double support_prob = 0.95)
    : values(Rcpp::as<values_t>(values))
    {
        init_bounds(support_prob);
    }

    ImportedValues(const Rcpp::MatrixColumn<REALSXP>& values, double support_prob = 0.95)
    : values(values.begin(), values.end())
    {
        init_bounds(support_prob);
    }

    ImportedValues(const values_t& values, double support_prob = 0.95)
    : values(values)
    {
        init_bounds(support_prob);
    }

    double defaultBinWidth(std::size_t nsteps) const {
        return (support_max - support_min) / nsteps;
    }
};
