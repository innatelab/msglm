#include <cmath>

// compresses the value,
// so that the result grows as sqrt(x) if x is above the threshold
// could be used to compress the p-values (-log10(p-value))
inline double sqrt_compress(double x,
                            double threshold = 10.0,
                            double hard_threshold_factor = 3)
{
    if (x < 0) throw std::invalid_argument("compressed value should be non-negative");
    if (x < threshold) return x;
    if (std::isfinite(x)) {
        const double t = x - threshold;
        return std::min(threshold + sqrt(t), threshold * hard_threshold_factor);
    } else {
        return threshold * hard_threshold_factor;
    }
};

inline double pvalue_sqrt_compress(double p,
                                   double mlog10_threshold = 10.0,
                                   double mlog10_hard_threshold_factor = 3)
{
    return exp10(-sqrt_compress(-log10(p), mlog10_threshold, mlog10_hard_threshold_factor));
};
