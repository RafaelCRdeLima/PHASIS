#include "phasis/trace.hpp"
#include "phasis/units.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace phasis {

namespace {

constexpr double kPiOver2 = 1.5707963267948966;

// Um dos dois ramos do raio: entrada (antes do ponto de retorno) e
// saida (depois). Para perfil esfericamente simetrico os dois dao o
// mesmo valor, mas NAO exploramos isso: na Fase 3 rho depende de theta
// e a simetria quebra. Cada ramo e integrado por conta propria.
enum class Branch { Incoming, Outgoing };

struct Integrand {
    const Metric&         metric;
    const DensityProfile& profile;
    const CrossSection&   xsec;

    double E_inf;
    double r_turn;
    double f_turn;
    double inclination;
    double psi_node;
    Branch branch;
    bool   spherical;

    // Angulo polar no ponto r.
    //
    // Em metrica esfericamente simetrica o movimento e PLANAR: basta
    // girar o plano orbital para cobrir qualquer geodesica. Um ponto a
    // angulo interno psi (medido do nodo ascendente) tem
    //
    //     cos(theta) = sin(inclinacao) * sin(psi)
    //
    // Para perfil esferico theta e irrelevante e devolvemos pi/2 sem
    // calcular psi -- que custaria uma quadratura aninhada por no.
    // O caminho nao esferico e da Fase 3.
    double theta_at(double r) const
    {
        if (spherical) return kPiOver2;
        (void)r;
        throw std::logic_error(
            "trace_ray: perfil nao esferico exige o angulo psi ao longo do "
            "raio, que e trabalho da Fase 3");
    }

    // Integrando na variavel s, com  r = r_turn + s^2.
    //
    //     dl/dr * dr/ds = sqrt(h(r)) / sqrt(w) * 2s
    //                   = 2 r sqrt( h(r) f(r_turn) / T(r) )
    //
    // (ver Metric::turning_factor: o s cancela analiticamente).
    double dl_ds(double r) const
    {
        const double T = metric.turning_factor(r, r_turn);
        if (!(T > 0.0)) return 0.0;
        return 2.0*r*std::sqrt(metric.h(r)*f_turn/T);
    }

    double rho_at(double r) const { return profile.rho(r, theta_at(r)); }

    // dtau/ds
    double operator()(double s) const
    {
        const double r = r_turn + s*s;

        const double rho = rho_at(r);
        if (rho <= 0.0) return 0.0;

        const double E_loc = metric.E_local(E_inf, r);
        const double sig   = xsec.sigma_tot(E_loc);

        return units::nucleons_per_gram*rho*sig*dl_ds(r);
    }
};

} // namespace

Result trace_ray(const Ray& ray,
                 const Metric& metric,
                 const DensityProfile& profile,
                 const CrossSection& xsec,
                 const IntegratorOpts& opts)
{
    Result out;

    const double b = ray.b_cm;
    if (!(b >= 0.0)) throw std::domain_error("trace_ray: b < 0");
    if (!(ray.E_inf_GeV > 0.0)) throw std::domain_error("trace_ray: E_inf <= 0");

    // ---- ponto de retorno --------------------------------------------
    // r_search_max = b: para f <= 1 vale r/sqrt(f) >= r, logo a raiz de
    // b = r/sqrt(f(r)) satisfaz r_turn <= b.
    const double r_turn = metric.r_turning(b, b);
    out.r_min_cm = r_turn;

    // ---- faixa util --------------------------------------------------
    // Cortar a faixa vazia aqui, em vez de deixar a quadratura descobrir
    // sozinha o degrau de rho, e o que mantem o custo baixo e a precisao
    // alta.
    const double r_lo = std::max(r_turn, profile.r_support_min());
    const double r_hi = profile.r_support_max();

    if (!(r_hi > r_lo)) {
        // O raio passa por fora da materia. tau = 0 exatamente, sem NaN.
        out.tau = 0.0;
        out.P_surv = 1.0;
        out.crosses_matter = false;
        return out;
    }
    out.crosses_matter = true;

    const double f_turn = metric.f(r_turn);
    if (!(f_turn > 0.0)) {
        throw std::domain_error("trace_ray: f(r_turn) <= 0");
    }

    const double s_lo = std::sqrt(std::max(0.0, r_lo - r_turn));
    const double s_hi = std::sqrt(r_hi - r_turn);

    // ---- os dois ramos -----------------------------------------------
    const Branch ramos[2] = { Branch::Incoming, Branch::Outgoing };

    for (Branch br : ramos) {

        Integrand ig{metric, profile, xsec,
                     ray.E_inf_GeV, r_turn, f_turn,
                     ray.inclination_rad, ray.psi_node_rad,
                     br, profile.is_spherical()};

        const QuadResult q = adaptive_simpson(
            [&](double s) { return ig(s); }, s_lo, s_hi, opts);

        out.tau     += q.value;
        out.n_evals += q.n_evals;
        if (q.depth_exhausted) out.tolerance_met = false;

        // diagnosticos: coluna de materia e comprimento proprio
        const QuadResult qx = adaptive_simpson(
            [&](double s) {
                const double r = r_turn + s*s;
                const double rho = ig.rho_at(r);
                return (rho > 0.0) ? rho*ig.dl_ds(r) : 0.0;
            }, s_lo, s_hi, opts);
        out.column_density += qx.value;

        const QuadResult ql = adaptive_simpson(
            [&](double s) {
                const double r = r_turn + s*s;
                return (ig.rho_at(r) > 0.0) ? ig.dl_ds(r) : 0.0;
            }, s_lo, s_hi, opts);
        out.path_length_cm += ql.value;
    }

    out.P_surv = std::exp(-out.tau);
    return out;
}

} // namespace phasis
