#include "phasis/trace.hpp"
#include "phasis/units.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace phasis {

namespace {

constexpr double kPiOver2 = 1.5707963267948966;
constexpr double kPi      = 3.141592653589793;
constexpr double kTwoPi   = 6.283185307179586;

// Ramos do raio. Em perfil esferico dao o mesmo valor, mas nao
// exploramos isso: na Fase 3 rho depende de theta e a simetria quebra.
enum class Branch { Incoming, Outgoing };

// =====================================================================
// Quadratura com tratamento do regime quase-critico.
//
// O integrando, apos r = r_t + s^2, vai como
//
//     g(s) ~ C / sqrt( T(r_t) + T'(r_t) s^2 )
//
// ou seja: constante para s << s_*, e C/(s sqrt(T')) para s >> s_*, com
//
//     s_* = sqrt( T(r_t) / T'(r_t) ).
//
// A cauda 1/s integrada de s_* a s_hi vale ln(s_hi/s_*) ~ -(1/2) ln T(r_t).
// Quando b -> b_crit+ tem-se T(r_t) -> 0 e isso DIVERGE logaritmicamente.
// A divergencia e FISICA: o raio enrola muitas vezes na esfera de fotons
// e o comprimento proprio realmente cresce sem limite. A substituicao
// s^2 matou a raiz quadrada, nao o logaritmo.
//
// Numericamente o problema e outro: o pico tem largura s_*, que pode ser
// 1e-10 de s_hi, e Simpson adaptativo em s gastaria dezenas de niveis
// para encontra-lo. Split em s_split = 8 s_* e, acima, substituicao
// exponencial s = s_split e^u -- onde ds/s = du e a cauda 1/s vira
// constante em u.
// =====================================================================
QuadResult integrate_peaked(const std::function<double(double)>& g,
                            double s_lo, double s_hi, double s_star,
                            const IntegratorOpts& opts)
{
    QuadResult out;
    if (!(s_hi > s_lo)) return out;

    const double s_split = 8.0*s_star;

    const bool pico_estreito = (s_star > 0.0) && (s_split < 0.1*s_hi);
    if (!pico_estreito) {
        return adaptive_simpson(g, s_lo, s_hi, opts);
    }

    auto acumula = [&out](const QuadResult& q) {
        out.value          += q.value;
        out.n_evals        += q.n_evals;
        out.max_depth_used  = std::max(out.max_depth_used, q.max_depth_used);
        out.depth_exhausted  = out.depth_exhausted  || q.depth_exhausted;
        out.budget_exhausted = out.budget_exhausted || q.budget_exhausted;
    };

    // (a) regiao do pico, direto em s
    const double hi1 = std::min(s_split, s_hi);
    if (hi1 > s_lo) acumula(adaptive_simpson(g, s_lo, hi1, opts));

    // (b) cauda, em u = ln(s/base)
    const double base = std::max(s_lo, s_split);
    if (s_hi > base) {
        auto gu = [&](double u) {
            const double s = base*std::exp(u);
            return g(s)*s;
        };
        acumula(adaptive_simpson(gu, 0.0, std::log(s_hi/base), opts));
    }

    return out;
}

// =====================================================================
struct Geometry {
    const Metric& metric;
    double r_turn;
    double f_turn;
    double T_turn;
    double T_slope;
    double s_star;

    // dl/dr * dr/ds = 2 r sqrt( h f(r_t) / T(r) ).  Ver Metric::turning_factor.
    double dl_ds(double r) const
    {
        const double T = metric.turning_factor(r, r_turn);
        if (!(T > 0.0)) return 0.0;
        return 2.0*r*std::sqrt(metric.h(r)*f_turn/T);
    }
};

struct Integrand {
    const Geometry&       geo;
    const DensityProfile& profile;
    const CrossSection&   xsec;
    double E_inf;
    Branch branch;
    bool   spherical;

    // Em metrica esfericamente simetrica o movimento e PLANAR, entao
    // basta girar o plano orbital. Para perfil esferico theta e
    // irrelevante. O caminho nao esferico e da Fase 3.
    double theta_at(double) const
    {
        if (spherical) return kPiOver2;
        throw std::logic_error(
            "trace_ray: perfil nao esferico exige psi ao longo do raio (Fase 3)");
    }

    double rho_at(double r) const { return profile.rho(r, theta_at(r)); }

    double operator()(double s) const
    {
        const double r = geo.r_turn + s*s;
        const double rho = rho_at(r);
        if (rho <= 0.0) return 0.0;

        // sigma DENTRO da integral, avaliada na energia LOCAL. Na Fase 1
        // f = 1 e isto fatoraria; nao fatoramos de proposito -- e o que
        // o teste T10 verifica.
        const double E_loc = geo.metric.E_local(E_inf, r);
        return units::nucleons_per_gram*rho*xsec.sigma_tot(E_loc)*geo.dl_ds(r);
    }
};

// =====================================================================
// Angulo azimutal total varrido, INT dphi/dr dr de r_t ate o infinito.
//
//   dphi/dr = (b/r^2) sqrt(h) / sqrt(1/f - b^2/r^2)
//           = (b/r^2) sqrt(h f) / sqrt(w)
//
// O sqrt(f) NAO pode ser esquecido: como 1/sqrt(f) ~ 1 + r_s/2r, deixa-lo
// de fora injeta um termo de ordem r_s/b em Delta_phi -- exatamente a
// ordem da propria resposta, dando 3 r_s/b em vez de 2 r_s/b.
// (Em Schwarzschild h f = 1 e o produto some; em metricas gerais nao.)
//
// Com u = 1/r e depois u = u_t - v^2 (u_t = 1/r_t) a integral fica
//
//   INT_0^{sqrt(u_t)} 2 b sqrt( h(r) f(r) u_t f(r_t) / ( u T(r) ) ) dv,  r = 1/u
//
// finita, sem cauda infinita e sem singularidade. Confirmacao analitica:
// em Minkowski isto da exatamente pi/2, logo Delta_phi = 2*(pi/2) - pi = 0.
// =====================================================================
double half_sweep_angle(const Geometry& geo, double b, const IntegratorOpts& opts)
{
    if (!(geo.r_turn > 0.0) || !(b > 0.0)) return kPiOver2;

    const double u_t = 1.0/geo.r_turn;

    auto integrando = [&](double v) {
        const double u = u_t - v*v;
        if (!(u > 0.0)) return 0.0;
        const double r = 1.0/u;
        const double T = geo.metric.turning_factor(r, geo.r_turn);
        if (!(T > 0.0)) return 0.0;
        return 2.0*b*std::sqrt(geo.metric.h(r)*geo.metric.f(r)
                               *u_t*geo.f_turn/(u*T));
    };

    // O pico quase-critico reaparece aqui com largura v_* = s_*/r_t:
    // perto de v = 0 vale T ~ T(r_t) + T'(r_t) r_t^2 v^2.
    const double v_star = geo.s_star/geo.r_turn;

    return integrate_peaked(integrando, 0.0, std::sqrt(u_t), v_star, opts).value;
}

} // namespace

// =====================================================================
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

    // ---- captura ------------------------------------------------------
    // Verificada ANTES de qualquer coisa, e por criterio proprio da
    // metrica -- nunca por "o root-finder nao convergiu".
    if (b > 0.0 && metric.is_captured(b)) {
        out.captured = true;
        out.tau      = std::numeric_limits<double>::infinity();
        out.P_surv   = 0.0;
        out.r_min_cm = 0.0;
        out.E_loc_max_GeV = std::numeric_limits<double>::infinity();
        return out;
    }

    // ---- ponto de retorno ---------------------------------------------
    const double r_turn = metric.r_turning(b, b);
    out.r_min_cm = r_turn;

    const double f_turn = metric.f(r_turn);
    if (!(f_turn > 0.0)) throw std::domain_error("trace_ray: f(r_turn) <= 0");

    const double T_turn  = metric.turning_factor_at_turn(r_turn);
    const double T_slope = metric.turning_factor_slope(r_turn);

    Geometry geo{metric, r_turn, f_turn, T_turn, T_slope, 0.0};
    geo.s_star = (T_turn > 0.0 && T_slope > 0.0)
               ? std::sqrt(T_turn/T_slope) : 0.0;

    // eps = T(r_t)/r_t, adimensional, mede a distancia a criticalidade.
    // So faz sentido com r_t > 0: em b = 0 tem-se T(r_t) = 0 sem que haja
    // qualquer singularidade (o fator 2r do integrando tambem se anula).
    if (r_turn > 0.0) {
        const double eps = T_turn/r_turn;
        out.near_critical = (eps <= opts.near_critical_eps);
    }

    // ---- deflexao ------------------------------------------------------
    // Calculada quando pedida, e sempre no regime quase-critico, onde o
    // numero de voltas e a informacao que interessa.
    if (opts.want_deflection || out.near_critical) {
        const double meia = half_sweep_angle(geo, b, opts);
        out.deflection_rad = 2.0*meia - kPi;
        // Voltas ALEM de uma linha reta: em Minkowski o angulo total
        // varrido de -inf a +inf e exatamente pi, entao a referencia e
        // Delta_phi = 0.
        out.winding_turns  = out.deflection_rad/kTwoPi;
    }

    // ---- faixa util ----------------------------------------------------
    const double r_lo = std::max(r_turn, profile.r_support_min());
    const double r_hi = profile.r_support_max();

    out.E_loc_max_GeV = metric.E_local(ray.E_inf_GeV,
                                       (r_hi > r_lo) ? r_lo : r_turn);

    if (!(r_hi > r_lo)) {
        out.tau = 0.0;
        out.P_surv = 1.0;
        out.crosses_matter = false;
        return out;
    }
    out.crosses_matter = true;

    const double s_lo = std::sqrt(std::max(0.0, r_lo - r_turn));
    const double s_hi = std::sqrt(r_hi - r_turn);

    // ---- os dois ramos --------------------------------------------------
    const Branch ramos[2] = { Branch::Incoming, Branch::Outgoing };

    for (Branch br : ramos) {
        Integrand ig{geo, profile, xsec, ray.E_inf_GeV, br, profile.is_spherical()};

        const QuadResult q =
            integrate_peaked([&](double s){ return ig(s); }, s_lo, s_hi, geo.s_star, opts);
        out.tau     += q.value;
        out.n_evals += q.n_evals;
        if (q.depth_exhausted || q.budget_exhausted) out.tolerance_met = false;

        const QuadResult qx = integrate_peaked([&](double s){
            const double r = geo.r_turn + s*s;
            const double rho = ig.rho_at(r);
            return (rho > 0.0) ? rho*geo.dl_ds(r) : 0.0;
        }, s_lo, s_hi, geo.s_star, opts);
        out.column_density += qx.value;

        const QuadResult ql = integrate_peaked([&](double s){
            const double r = geo.r_turn + s*s;
            return (ig.rho_at(r) > 0.0) ? geo.dl_ds(r) : 0.0;
        }, s_lo, s_hi, geo.s_star, opts);
        out.path_length_cm += ql.value;
    }

    out.P_surv = std::exp(-out.tau);
    return out;
}

} // namespace phasis
