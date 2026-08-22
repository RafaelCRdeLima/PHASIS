#include "phasis/density.hpp"

#include <cmath>

namespace phasis {

double PowerLawHalo::rho(double r, double) const
{
    if (r < r_in_ || r > r_out_) return 0.0;
    return rho0_*std::pow(r/r0_, -p_);
}

double FlaredThinDisk::H_of_R(double R) const
{
    return H0_*std::pow(R/R0_, q_);
}

double FlaredThinDisk::rho(double r, double theta) const
{
    if (r < r_in_ || r > r_out_) return 0.0;

    const double R = r*std::sin(theta);
    const double z = r*std::cos(theta);

    if (!(R > 0.0)) return 0.0;    // eixo: R = 0, sem materia

    const double H = H_of_R(R);
    if (!(H > 0.0)) return 0.0;

    // Em log, e nao direto. Perto do eixo (R/R_0)^(-p) explode enquanto a
    // gaussiana colapsa; o produto tende a zero, mas calculado
    // separadamente da inf*0. Em log a soma e bem comportada e o corte
    // em -700 devolve zero antes de qualquer overflow.
    const double ln_rho = std::log(rho0_)
                        - p_*std::log(R/R0_)
                        - 0.5*(z/H)*(z/H);

    if (ln_rho < -700.0) return 0.0;
    return std::exp(ln_rho);
}

double QuasiSphericalADAF::rho(double r, double) const
{
    if (r < r_in_ || r > r_out_) return 0.0;
    return rho0_*std::pow(r/R0_, -1.5);
}

} // namespace phasis
