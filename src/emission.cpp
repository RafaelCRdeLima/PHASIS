#include "phasis/emission.hpp"

#include <cmath>
#include <stdexcept>

namespace phasis {

namespace {
constexpr double kPiOver2 = 1.5707963267948966;
}

double b_from_emission(double r_emit, double psi, const Metric& metric)
{
    if (!(r_emit > 0.0)) throw std::runtime_error("b_from_emission: r_emit <= 0");
    const double f = metric.f(r_emit);
    if (!(f > 0.0)) throw std::runtime_error("b_from_emission: r_emit dentro do horizonte");
    return r_emit*std::fabs(std::sin(psi))/std::sqrt(f);
}

const char* to_string(Topology t)
{
    switch (t) {
        case Topology::OutboundOnly: return "OutboundOnly";
        case Topology::Turning:      return "Turning";
        case Topology::Captured:     return "Captured";
    }
    return "?";
}

Topology classify(double r_emit, double psi, const Metric& metric, double* b_out)
{
    const double b = b_from_emission(r_emit, psi, metric);
    if (b_out) *b_out = b;

    const bool para_fora = (psi < kPiOver2);
    const double r_ph    = metric.r_photon_sphere();

    // has_turning: existe raiz externa de g(r) = b, ou seja b > b_crit.
    // Em Minkowski b_crit = 0 e isto e sempre verdade para b > 0.
    const bool tem_retorno = (b > 0.0) && !metric.is_captured(b);

    if (r_emit > r_ph) {
        if (para_fora) return Topology::OutboundOnly;      // sempre escapa
        return tem_retorno ? Topology::Turning : Topology::Captured;
    }

    // Abaixo da esfera de fotons g DECRESCE: a condicao inverte.
    if (para_fora) return tem_retorno ? Topology::Captured : Topology::OutboundOnly;
    return Topology::Captured;
}

double escape_cone_angle(double r_emit, const Metric& metric)
{
    const double r_ph = metric.r_photon_sphere();
    if (!(r_ph > 0.0)) return kPiOver2;      // sem esfera de fotons: nada captura

    // b_crit = g(r_ph) = r_ph/sqrt(f(r_ph))
    const double b_crit = r_ph/std::sqrt(metric.f(r_ph));

    double s = b_crit*std::sqrt(metric.f(r_emit))/r_emit;
    if (s > 1.0) s = 1.0;
    if (s < 0.0) s = 0.0;
    return std::asin(s);
}

double captured_fraction(double r_emit, const Metric& metric)
{
    const double r_ph = metric.r_photon_sphere();
    if (!(r_ph > 0.0)) return 0.0;
    if (r_emit <= r_ph) {
        // Dentro da esfera de fotons: captura tudo que vai para dentro,
        // mais o que vai para fora com b > b_crit. Em r = r_ph exato isso
        // da exatamente metade do ceu.
        const double psi_c = escape_cone_angle(r_emit, metric);
        return 1.0 - (1.0 - std::cos(psi_c))/2.0;
    }
    return (1.0 - std::cos(escape_cone_angle(r_emit, metric)))/2.0;
}

Ray emit_ray(double r_emit, double psi, double E_inf_GeV,
             const Metric& metric, double inclination_rad, double psi_turn_rad)
{
    Ray r;
    r.E_inf_GeV       = E_inf_GeV;
    r.b_cm            = b_from_emission(r_emit, psi, metric);
    r.inclination_rad = inclination_rad;
    r.psi_turn_rad    = psi_turn_rad;
    r.r_emit_cm       = r_emit;
    r.outward         = (psi < kPiOver2);
    return r;
}

} // namespace phasis
