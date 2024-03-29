#include "imported_values.h"

#include "logging.h"

#include <boost/array.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/extended_p_square.hpp>

namespace bacc = boost::accumulators;

void ImportedValues::init_bounds(double support_prob)
{
    if (support_prob < 0.5 || support_prob > 1.0) THROW_EXCEPTION(std::invalid_argument,
            "support_prob=%g, should be in [0.5, 1.0] range", support_prob);
    double valmin = !values.empty() ? values[0] : std::numeric_limits<value_t>::quiet_NaN();
    double valmax = !values.empty() ? values[0] : std::numeric_limits<value_t>::quiet_NaN();
    const double tail_prob = 0.5*(1.0 - support_prob);
    bacc::accumulator_set<double, bacc::stats<bacc::tag::extended_p_square> >
        quant_acc(bacc::tag::extended_p_square::probabilities = boost::array<double, 2>{tail_prob, 1.0 - tail_prob});
    if (!values.empty()) quant_acc(values[0]);
    for ( std::size_t i = 1; i < values.size(); ++i ) {
      const double val = values[i];
      quant_acc(val);
      if (valmin > val) valmin = val;
      if (valmax < val) valmax = val;
    }
    val_min = valmin;
    val_max = valmax;
    support_min = bacc::extended_p_square(quant_acc)[0];
    support_max = bacc::extended_p_square(quant_acc)[1];
    LOG_DEBUG2("minmax=[%g, %g] support=[%g, %g]", val_min, val_max, support_min, support_max);
}
