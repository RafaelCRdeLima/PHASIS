#include "phasis/density.hpp"

#include <cmath>

namespace phasis {

double PowerLawHalo::rho(double r, double) const
{
    if (r < r_in_ || r > r_out_) return 0.0;
    return rho0_*std::pow(r/r0_, -p_);
}

} // namespace phasis
