#include "phasis/trace.hpp"
#include "phasis/emission.hpp"
#include "phasis/ode.hpp"
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

    // Modo da substituicao r = r_anchor + s^2.
    //
    //  Normal : ha ponto de retorno e o raio o alcanca. r_anchor = r_t,
    //           e a raiz de w em r_t e removida por turning_factor.
    //  NoTurn : NAO ha ponto de retorno na faixa percorrida -- e o caso
    //           de um raio emitido PARA FORA com b < b_crit, que escapa.
    //           Aqui w = 1 - f b^2/r^2 nunca se anula, entao nao ha
    //           singularidade e nao ha o que fatorar: calcula-se direto.
    //           r_anchor = r_emit.
    //  Radial : b = 0. dl/dr = sqrt(h) exatamente; r_t e f(r_t) nem sao
    //           invocados (em Schwarzschild seriam singulares).
    enum class Mode { Normal, NoTurn, Radial };
    Mode   mode;
    double r_anchor;
    double b;

    // dl/dr * dr/ds = 2 r sqrt( h f(r_t) / T(r) ).  Ver Metric::turning_factor.
    //
    // O caso RADIAL (b = 0) vai por caminho separado: la r_t = 0 e a
    // fatoracao por turning_factor fica mal-condicionada -- em
    // Schwarzschild ela nem esta definida em r_t = 0. Radialmente
    // dl/dr = sqrt(h) exatamente, e com r = s^2 isso da
    //     dl/ds = 2 sqrt(r) sqrt(h(r)),
    // que e o limite continuo da formula geral (checado em T34).
    double s_de_r(double r) const { return std::sqrt(std::max(0.0, r - r_anchor)); }
    double r_de_s(double s) const { return r_anchor + s*s; }

    double dl_ds(double r) const
    {
        const double h = metric.h(r);
        const double s = s_de_r(r);

        if (mode == Mode::Radial) return 2.0*s*std::sqrt(h);

        if (mode == Mode::NoTurn) {
            const double w = 1.0 - metric.f(r)*b*b/(r*r);
            if (!(w > 0.0)) return 0.0;
            return 2.0*s*std::sqrt(h/w);
        }

        const double T = metric.turning_factor(r, r_turn);
        if (!(T > 0.0)) return 0.0;
        return 2.0*r*std::sqrt(h*f_turn/T);
    }
};

struct Integrand {
    const Geometry&       geo;
    const DensityProfile& profile;
    const CrossSection&   xsec;
    double E_inf;
    Branch branch;
    bool   spherical;
    bool   freeze;   // sigma avaliada em E_inf, nao em E_loc

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
        const double r = geo.r_de_s(s);
        const double rho = rho_at(r);
        if (rho <= 0.0) return 0.0;

        // sigma DENTRO da integral, avaliada na energia LOCAL. Na Fase 1
        // f = 1 e isto fatoraria; nao fatoramos de proposito -- e o que
        // o teste T10 verifica.
        const double E_loc = freeze ? E_inf : geo.metric.E_local(E_inf, r);
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


RayDerivatives ray_derivatives(const Metric& metric,
                               double r, double r_turn, double b)
{
    RayDerivatives d;

    const double h = metric.h(r);

    if (b == 0.0) {
        // Radial: dl/dr = sqrt(h), e dpsi/ds = 0 por definicao.
        d.dl_ds   = 2.0*std::sqrt(r)*std::sqrt(h);
        d.dpsi_ds = 0.0;
        return d;
    }

    const double T = metric.turning_factor(r, r_turn);
    if (!(T > 0.0)) return d;

    const double fr  = metric.f(r);
    const double f_t = metric.f(r_turn);

    d.dl_ds   = 2.0*r*std::sqrt(h*f_t/T);
    d.dpsi_ds = (2.0*b/r)*std::sqrt(h*fr*f_t/T);
    return d;
}

// =====================================================================
// Rota de EDO: obrigatoria quando rho depende de theta.
//
// O acoplamento: psi(r) e necessario DENTRO do integrando de tau, entao
// as duas quadraturas deixam de ser independentes. Resolver com
// sub-quadratura aninhada custaria O(N^2) por raio -- inviavel para 1e5
// raios. Aqui o estado avanca junto, numa passada so.
//
// Com a mesma fatoracao das Fases 1-2, sqrt(w) = s sqrt(T)/(r sqrt(f_t))
// e r = r_t + s^2, os dois integrandos ficam regulares em s:
//
//     dl/ds   = 2 r     sqrt( h   f_t / T )
//     dpsi/ds = (2 b/r) sqrt( h f f_t / T )
//
// Identidade util, verificada em T13:  dpsi/dl = b sqrt(f) / r^2.
//
// Estado:
//   y[0] psi_off   deslocamento de psi a partir do ponto de retorno (>= 0)
//   y[1] l         comprimento proprio de UM ramo
//   y[2] tau_in    ramo de entrada, psi = psi_t - psi_off
//   y[3] tau_out   ramo de saida,   psi = psi_t + psi_off
//   y[4] X_in      coluna, g/cm^2
//   y[5] X_out
//
// Os dois ramos avancam no MESMO vetor de estado porque compartilham
// psi_off e l; so o sinal com que psi_off entra em theta difere. Nunca
// se dobra um ramo -- T15 existe para pegar exatamente esse bug.
// =====================================================================
namespace {

struct OdeCtx {
    const Geometry&       geo;
    const DensityProfile& profile;
    const CrossSection&   xsec;
    double E_inf;
    double b;
    double sin_i;
    double psi_t;
    bool   freeze;

    // Trecho SEM materia: so psi e l avancam.
    //
    // Nao e otimizacao. Os perfis tem corte duro em r_in, e o ultimo
    // estagio de Runge-Kutta do trecho [0, s_lo] cai exatamente em
    // r = r_in, onde rho salta de 0 para finito. Um degrau bem na
    // fronteira faz o estimador de erro explodir, o passo encolher ate
    // h_min, e encolher nao ajuda -- a descontinuidade esta no extremo.
    // Zerando as derivadas de tau e X onde nao ha materia por
    // construcao, o degrau some do controle de erro.
    bool matter = true;

    double theta_de(double psi) const {
        double c = sin_i*std::sin(psi);
        if (c >  1.0) c =  1.0;
        if (c < -1.0) c = -1.0;
        return std::acos(c);
    }

    void operator()(double s, const OdeState<6>& y, OdeState<6>& dy) const
    {
        const double r = geo.r_turn + s*s;

        const RayDerivatives d = ray_derivatives(geo.metric, r, geo.r_turn, b);
        if (!(d.dl_ds > 0.0)) { dy.fill(0.0); return; }

        const double dl   = d.dl_ds;
        const double dpsi = d.dpsi_ds;

        const double E_loc = matter
            ? (freeze ? E_inf : geo.metric.E_local(E_inf, r))
            : 0.0;
        const double sig   = matter ? xsec.sigma_tot(E_loc) : 0.0;

        dy[0] = dpsi;
        dy[1] = dl;

        if (!matter) { dy[2] = dy[3] = dy[4] = dy[5] = 0.0; return; }

        const double rho_in  = profile.rho(r, theta_de(psi_t - y[0]));
        const double rho_out = profile.rho(r, theta_de(psi_t + y[0]));

        dy[2] = units::nucleons_per_gram*rho_in *sig*dl;
        dy[3] = units::nucleons_per_gram*rho_out*sig*dl;
        dy[4] = rho_in *dl;
        dy[5] = rho_out*dl;
    }
};

} // namespace

static void trace_via_ode(Result& out,
                          const Ray& ray,
                          const Geometry& geo,
                          const DensityProfile& profile,
                          const CrossSection& xsec,
                          double s_lo,
                          double s_hi,
                          const IntegratorOpts& opts,
                          bool freeze)
{
    OdeCtx ctx{geo, profile, xsec, ray.E_inf_GeV, ray.b_cm,
               std::sin(ray.inclination_rad), ray.psi_turn_rad, freeze};

    OdeState<6> y{};   // tudo zero: psi_off, l, tau, X partem do ponto de retorno

    OdeOpts oo;
    oo.rel_tol   = opts.ode_rel_tol;
    oo.max_steps = opts.ode_max_steps;

    auto passo = [&](double a, double b, bool com_materia) {
        if (!(b > a)) return;
        ctx.matter = com_materia;
        const OdeStats st = dopri54<6>(
            [&](double s, const OdeState<6>& yy, OdeState<6>& dd){ ctx(s, yy, dd); },
            a, b, y, oo);
        if (!st.ok) out.tolerance_met = false;
        out.n_evals += st.n_evals;
    };

    // Comeca em s = 0 SEMPRE, mesmo que a materia so comece mais acima:
    // psi precisa acumular desde o ponto de retorno, senao o mapeamento
    // psi -> theta fica deslocado.
    //
    // MAS o intervalo e CORTADO na fronteira interna do suporte. Sem
    // isso, no trecho vazio as derivadas de tau sao identicamente nulas,
    // o controlador so ve psi e l (que sao suaves), o passo cresce 5x por
    // iteracao e o solver salta por cima da materia. Foi exatamente o que
    // aconteceu em T12 com b/r_s = 3, onde r_t fica ABAIXO de r_in:
    // tau saiu zero. O corte forca uma fronteira de passo onde a materia
    // comeca.
    passo(0.0,  s_lo, false);   // vazio por construcao
    passo(s_lo, s_hi, true);

    out.tau_inbound     = y[2];
    out.tau_outbound    = y[3];
    out.column_inbound  = y[4];
    out.column_outbound = y[5];

    out.tau            = y[2] + y[3];
    out.column_density = y[4] + y[5];
    out.path_length_cm = 2.0*y[1];

    // ---- extremos de theta e cruzamentos do plano, analiticamente ------
    // psi percorre [psi_t - psi_max, psi_t + psi_max], e
    // cos(theta) = sin(i) sin(psi).
    const double psi_max = y[0];
    const double lo = ctx.psi_t - psi_max;
    const double hi = ctx.psi_t + psi_max;

    double smin = std::min(std::sin(lo), std::sin(hi));
    double smax = std::max(std::sin(lo), std::sin(hi));
    // extremos internos de sin: psi = pi/2 + k pi
    for (double k = std::floor((lo - kPiOver2)/kPi); k <= std::ceil((hi - kPiOver2)/kPi); ++k) {
        const double psi = kPiOver2 + k*kPi;
        if (psi >= lo && psi <= hi) {
            smin = std::min(smin, std::sin(psi));
            smax = std::max(smax, std::sin(psi));
        }
    }
    double cmin = ctx.sin_i*smin, cmax = ctx.sin_i*smax;
    if (cmin > cmax) std::swap(cmin, cmax);
    cmin = std::max(-1.0, cmin);
    cmax = std::min( 1.0, cmax);
    out.theta_min_rad = std::acos(cmax);
    out.theta_max_rad = std::acos(cmin);

    // z = 0  <=>  cos(theta) = 0  <=>  sin(psi) = 0  <=>  psi = k pi.
    // Com i = 0 o raio ESTA no plano o tempo todo: nao cruza.
    if (std::fabs(ctx.sin_i) > 1.0e-15) {
        const long k0 = static_cast<long>(std::ceil (lo/kPi));
        const long k1 = static_cast<long>(std::floor(hi/kPi));
        out.n_disk_crossings = static_cast<int>(std::max(0L, k1 - k0 + 1));
    }
}

// =====================================================================
Result trace_ray(const Ray& ray,
                 const Metric& metric,
                 const DensityProfile& profile,
                 const CrossSection& xsec,
                 const IntegratorOpts& opts,
                 bool freeze_redshift)
{
    Result out;

    const double b = ray.b_cm;
    if (!(b >= 0.0)) throw std::domain_error("trace_ray: b < 0");
    if (!(ray.E_inf_GeV > 0.0)) throw std::domain_error("trace_ray: E_inf <= 0");

    // ---- topologia ----------------------------------------------------
    // Para raio EMITIDO a distancia finita, o criterio de escape depende
    // de que lado da esfera de fotons r_emit esta -- sao quatro casos,
    // nao tres. Ver emission.hpp.
    //
    // Metric::is_captured significa "nao existe ponto de retorno", que
    // coincide com captura SO para raio vindo do infinito.
    const bool emitido = (ray.r_emit_cm > 0.0);
    Topology topo = Topology::Turning;
    if (emitido) {
        const double f_e = metric.f(ray.r_emit_cm);
        const double sin_psi = (f_e > 0.0)
            ? std::min(1.0, b*std::sqrt(f_e)/ray.r_emit_cm) : 0.0;
        const double psi = ray.outward ? std::asin(sin_psi) : (kPi - std::asin(sin_psi));
        topo = classify(ray.r_emit_cm, psi, metric);
        if (topo == Topology::Captured) {
            out.captured = true;
            out.tau      = std::numeric_limits<double>::infinity();
            out.P_surv   = 0.0;
            out.E_loc_max_GeV = std::numeric_limits<double>::infinity();
            return out;
        }
    }

    // ---- captura (raio vindo do infinito) -----------------------------
    //
    // SEM o guarda `b > 0`: um raio RADIAL vindo do infinito cai no
    // buraco negro, e Schwarzschild::is_captured(0) devolve true
    // corretamente. Com o guarda, b = 0 pulava a checagem e ia direto
    // para r_turning, que entao lancava. Minkowski::is_captured(0)
    // devolve false, entao o caso plano segue funcionando.
    if (!emitido && metric.is_captured(b)) {
        out.captured = true;
        out.tau      = std::numeric_limits<double>::infinity();
        out.P_surv   = 0.0;
        out.r_min_cm = 0.0;
        out.E_loc_max_GeV = std::numeric_limits<double>::infinity();
        return out;
    }

    // ---- ponto de retorno ---------------------------------------------
    //
    // Raio RADIAL (b = 0): nao ha ponto de retorno a invocar. Em
    // Schwarzschild r_turning(0) lancaria, e f(0) e singular -- mas o
    // caminho radial nao usa nenhum dos dois: dl/dr = sqrt(h) direto.
    // Fixamos r_turn = 0 e f_turn = 1 como valores inertes.
    const bool radial = (b == 0.0);

    const bool sem_retorno = radial || metric.is_captured(b);
    const double r_turn = sem_retorno ? 0.0 : metric.r_turning(b, b);
    out.r_min_cm = sem_retorno ? (emitido ? ray.r_emit_cm : 0.0) : r_turn;

    const double f_turn = sem_retorno ? 1.0 : metric.f(r_turn);
    if (!sem_retorno && !(f_turn > 0.0))
        throw std::domain_error("trace_ray: f(r_turn) <= 0");

    const double T_turn  = sem_retorno ? 0.0 : metric.turning_factor_at_turn(r_turn);
    const double T_slope = sem_retorno ? 1.0 : metric.turning_factor_slope(r_turn);

    // Escolha do modo. NoTurn e o caso que a especificacao original nao
    // previa: um raio emitido PARA FORA com b < b_crit nao tem ponto de
    // retorno nenhum -- nao um ponto de retorno "virtual", mas ausencia.
    Geometry::Mode modo = Geometry::Mode::Normal;
    double r_anchor = r_turn;
    if (radial) {
        modo = Geometry::Mode::Radial;
        r_anchor = (emitido && ray.outward) ? ray.r_emit_cm : 0.0;
    } else if (emitido && topo == Topology::OutboundOnly && sem_retorno) {
        modo = Geometry::Mode::NoTurn;
        r_anchor = ray.r_emit_cm;
    }

    Geometry geo{metric, r_turn, f_turn, T_turn, T_slope, 0.0, modo, r_anchor, b};
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
    const double r_lo = std::max(r_anchor, profile.r_support_min());
    const double r_hi = profile.r_support_max();

    // Ponto mais fundo REALMENTE percorrido. Para um raio emitido para
    // fora isso e r_emit, nao o ponto de retorno virtual (que pode ser 0
    // e onde f seria singular).
    double r_fundo = (r_hi > r_lo) ? r_lo : r_turn;
    if (emitido && topo == Topology::OutboundOnly) {
        r_fundo = std::max(r_fundo, ray.r_emit_cm);
    }
    if (!(r_fundo > 0.0)) r_fundo = std::max(profile.r_support_min(), 1.0);
    out.E_loc_max_GeV = metric.E_local(ray.E_inf_GeV, r_fundo);

    if (!(r_hi > r_lo)) {
        out.tau = 0.0;
        out.P_surv = 1.0;
        out.crosses_matter = false;
        return out;
    }
    out.crosses_matter = true;

    const double s_lo = geo.s_de_r(r_lo);
    const double s_hi = geo.s_de_r(r_hi);

    // ---- despacho -------------------------------------------------------
    // Perfil esferico: quadratura das Fases 1-2, validada em 1e-16.
    // Perfil com theta: EDO acoplada, obrigatoria.
    // force_ode: T12 confronta as duas no mesmo caso esferico.
    if (!profile.is_spherical() || opts.force_ode) {
        trace_via_ode(out, ray, geo, profile, xsec, s_lo, s_hi, opts,
                      freeze_redshift);
        out.P_surv = std::exp(-out.tau);
        return out;
    }

    // ---- faixas POR RAMO ------------------------------------------------
    // Raio do infinito: os dois ramos cobrem [s_lo, s_hi].
    // Emitido, Turning: entrada so vai de r_t ate r_emit.
    // Emitido, OutboundOnly: um ramo so, de r_emit para fora. Note que
    //   r_t continua definido -- e o ponto de retorno VIRTUAL, que o raio
    //   nunca atinge, e a fatoracao T(r) segue valida.
    double s_in_hi = s_hi, s_out_lo = s_lo;
    bool tem_entrada = true;

    if (emitido) {
        const double s_emit = geo.s_de_r(ray.r_emit_cm);
        if (topo == Topology::OutboundOnly) {
            tem_entrada = false;
            s_out_lo = std::max(s_lo, s_emit);
        } else {
            s_in_hi = std::min(s_hi, s_emit);
        }
    }

    // ---- os dois ramos, por quadratura ----------------------------------
    const Branch ramos[2] = { Branch::Incoming, Branch::Outgoing };

    double tau_ramo[2] = {0.0, 0.0};
    double col_ramo[2] = {0.0, 0.0};

    for (int ib = 0; ib < 2; ++ib) {
        if (ib == 0 && !tem_entrada) continue;
        const double a = (ib == 0) ? s_lo     : s_out_lo;
        const double c = (ib == 0) ? s_in_hi  : s_hi;
        if (!(c > a)) continue;
        Integrand ig{geo, profile, xsec, ray.E_inf_GeV, ramos[ib],
                     profile.is_spherical(), freeze_redshift};

        const QuadResult q =
            integrate_peaked([&](double s){ return ig(s); }, a, c, geo.s_star, opts);
        tau_ramo[ib] = q.value;
        out.n_evals += q.n_evals;
        if (q.depth_exhausted || q.budget_exhausted) out.tolerance_met = false;

        const QuadResult qx = integrate_peaked([&](double s){
            const double r = geo.r_de_s(s);
            const double rho = ig.rho_at(r);
            return (rho > 0.0) ? rho*geo.dl_ds(r) : 0.0;
        }, a, c, geo.s_star, opts);
        col_ramo[ib] = qx.value;

        const QuadResult ql = integrate_peaked([&](double s){
            const double r = geo.r_de_s(s);
            return (ig.rho_at(r) > 0.0) ? geo.dl_ds(r) : 0.0;
        }, a, c, geo.s_star, opts);
        out.path_length_cm += ql.value;
    }

    out.tau_inbound     = tau_ramo[0];
    out.tau_outbound    = tau_ramo[1];
    out.column_inbound  = col_ramo[0];
    out.column_outbound = col_ramo[1];
    out.tau             = tau_ramo[0] + tau_ramo[1];
    out.column_density  = col_ramo[0] + col_ramo[1];

    // perfil esferico: theta nao e visitado, e nao ha "cruzamento"
    out.theta_min_rad = kPiOver2;
    out.theta_max_rad = kPiOver2;

    out.P_surv = std::exp(-out.tau);
    return out;
}

} // namespace phasis
