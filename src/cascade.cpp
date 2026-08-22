#include "phasis/cascade.hpp"
#include "phasis/integrate.hpp"
#include "phasis/ode.hpp"
#include "phasis/trace.hpp"
#include "phasis/units.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace phasis {

// =====================================================================
EnergyGrid::EnergyGrid(double E_min, double E_max, int n_bins)
    : E_min_(E_min), E_max_(E_max), n_(n_bins)
{
    if (!(E_min > 0.0) || !(E_max > E_min)) {
        throw std::runtime_error("EnergyGrid: faixa de energia invalida");
    }
    if (n_bins < 1) throw std::runtime_error("EnergyGrid: n_bins < 1");

    ln_min_ = std::log(E_min);
    dln_    = (std::log(E_max) - ln_min_)/n_bins;
    ratio_  = std::exp(dln_);
}

double EnergyGrid::edge(int i) const
{
    return std::exp(ln_min_ + i*dln_);
}

double EnergyGrid::center(int i) const
{
    return std::exp(ln_min_ + (i + 0.5)*dln_);
}

double EnergyGrid::powerlaw_bin_integral(int i, double gamma) const
{
    const double a = edge(i), b = edge(i + 1);
    if (std::fabs(gamma - 1.0) < 1.0e-14) return std::log(b/a);
    return (std::pow(b, 1.0 - gamma) - std::pow(a, 1.0 - gamma))/(1.0 - gamma);
}

// =====================================================================
namespace {

// Integral de dsigma/dy normalizada, sobre [y0, y1], por Simpson
// composto com um numero fixo (grande) de pontos.
//
// Nao e adaptativa de proposito: a identidade de consistencia exige que
// os intervalos ADJACENTES usem a MESMA regra, senao os erros nao se
// cancelam na soma e a identidade fecha pior do que a precisao de cada
// pedaco. Com regra fixa por intervalo, a soma dos pedacos e uma regra
// composta legitima sobre [0,1].
// Quadratura ADAPTATIVA por pedaco, nao regra fixa.
//
// A tentativa inicial foi Simpson composto fixo, com o argumento de que
// a soma dos pedacos vira uma regra composta legitima sobre [0,1] e os
// erros cancelam. Isso vale para integrando suave -- ConstantY fechava
// em 5.6e-14 -- mas falha feio num kernel picado: com g ~ y^-1 perto de
// y_min, Simpson LINEAR errava por 1.6e-4 e a identidade nao fechava.
//
// Com adaptativa e tolerancia ABSOLUTA por pedaco, cada integral
// converge por si; a soma de N pedacos acumula no maximo N*abs_tol, que
// para N ~ 1e3 e abs_tol ~ 1e-16 sigma fica em 1e-13 relativo. Nao
// depende de cancelamento entre vizinhos.
double simpson_liso(const DifferentialCrossSection& xs, double E,
                    double y0, double y1, double abs_tol)
{
    if (!(y1 > y0)) return 0.0;

    IntegratorOpts o;
    o.rel_tol   = 1.0e-14;
    o.abs_tol   = abs_tol;
    o.max_depth = 40;
    o.max_evals = 2000000;

    return adaptive_simpson([&](double y){ return xs.dsigma_nc_dy(E, y); },
                            y0, y1, o).value;
}

// Integral sobre [y0,y1] partida nos pontos de quebra declarados pelo
// kernel. Sem isso, uma descontinuidade DENTRO de um painel de Simpson
// impede a identidade de consistencia de fechar: o erro nao cancela
// entre intervalos adjacentes.
double integra_g(const DifferentialCrossSection& xs, double E,
                 double y0, double y1, double abs_tol,
                 const std::vector<double>& quebras)
{
    if (!(y1 > y0)) return 0.0;

    std::vector<double> pts{y0};
    for (double q : quebras) if (q > y0 && q < y1) pts.push_back(q);
    pts.push_back(y1);
    std::sort(pts.begin(), pts.end());

    double total = 0.0;
    for (std::size_t k = 0; k + 1 < pts.size(); ++k) {
        total += simpson_liso(xs, E, pts[k], pts[k+1], abs_tol);
    }
    return total;
}

} // namespace

// =====================================================================
PowerLawY::PowerLawY(double sigma_nc0, double sigma_cc0, double E0,
                     double alpha_nc, double alpha_cc,
                     double beta, double y_min)
    : nc0_(sigma_nc0), cc0_(sigma_cc0), E0_(E0),
      anc_(alpha_nc), acc_(alpha_cc), beta_(beta), y_min_(y_min)
{
    if (!(y_min > 0.0) || !(y_min < 1.0)) {
        throw std::runtime_error("PowerLawY: y_min tem de estar em (0,1)");
    }
    // norm = integral_{y_min}^1 y^-beta dy
    norm_ = (std::fabs(beta - 1.0) < 1.0e-14)
          ? -std::log(y_min)
          : (1.0 - std::pow(y_min, 1.0 - beta))/(1.0 - beta);
}

double PowerLawY::dsigma_nc_dy(double E, double y) const
{
    if (y < y_min_ || y > 1.0) return 0.0;
    return sigma_nc(E)*std::pow(y, -beta_)/norm_;
}

double PowerLawY::Z_analytic(double gamma, int n_quad) const
{
    // integral_{y_min}^1 y^-beta (1-y)^(gamma-1) dy / norm,
    // em log y, que resolve o pico em y_min.
    if (n_quad % 2 != 0) ++n_quad;
    const double a = std::log(y_min_), b = 0.0;
    const double h = (b - a)/n_quad;
    auto f = [&](double t) {
        const double y = std::exp(t);
        return std::pow(y, 1.0 - beta_)*std::pow(1.0 - y, gamma - 1.0);
    };
    double soma = f(a) + f(b);
    for (int k = 1; k < n_quad; ++k) soma += (k % 2 ? 4.0 : 2.0)*f(a + k*h);
    return (soma*h/3.0)/norm_;
}

// =====================================================================
CascadeKernel::CascadeKernel(const EnergyGrid& grid,
                             const DifferentialCrossSection& xsec,
                             double consistency_tol,
                             int /*n_y_quad, obsoleto: a quadratura e adaptativa*/)
    : grid_(grid), xsec_(xsec), n_(grid.n())
{
    if (!xsec.shape_is_E_independent()) {
        throw std::runtime_error(
            "CascadeKernel: esta fase so implementa kernels com forma de y "
            "independente de E (dsigma/dy = sigma_nc(E) g(y)). Um kernel com "
            "forma dependente de E exige K_ij(E) tabelada -- ver Fase 4b.");
    }

    G_.assign(static_cast<std::size_t>(n_)*static_cast<std::size_t>(n_), 0.0);
    G_leak_.assign(static_cast<std::size_t>(n_), 0.0);
    residual_.assign(static_cast<std::size_t>(n_), 0.0);

    // Energia de referencia: como a forma independe de E, qualquer uma
    // serve para extrair g(y). Normalizamos por sigma_nc na MESMA energia.
    const double E_ref = grid.center(0);
    const double s_ref = xsec.sigma_nc(E_ref);

    // Caso degenerado: sem NC nao ha o que redistribuir. A identidade
    // NORMALIZADA (soma G = 1) e 0/0 e nao faz sentido; a identidade
    // real, soma_i K_ij + leak_j = sigma_nc = 0, vale trivialmente.
    // G fica tudo zero e o residuo e zero, nao -1.
    if (!(s_ref > 0.0)) {
        no_regeneration_ = true;
        return;
    }
    const std::vector<double> quebras = xsec.y_breakpoints();

    // Tolerancia ABSOLUTA por pedaco, escalada por sigma_nc. E o que
    // garante que a soma de muitos pedacos ainda feche a identidade.
    const double abs_piece = 1.0e-16*std::fabs(s_ref);

    for (int j = 0; j < n_; ++j) {
        const double Ej = grid.center(j);

        double soma = 0.0;

        for (int i = 0; i <= j; ++i) {
            // E = (1-y) Ej cai no bin i  <=>  y em
            //   [ 1 - edge(i+1)/Ej , 1 - edge(i)/Ej ]
            double y0 = 1.0 - grid.edge(i + 1)/Ej;
            double y1 = 1.0 - grid.edge(i)/Ej;
            y0 = std::max(0.0, y0);
            y1 = std::min(1.0, y1);
            if (!(y1 > y0)) continue;

            const double v = (s_ref > 0.0)
                           ? integra_g(xsec, E_ref, y0, y1, abs_piece, quebras)/s_ref
                           : 0.0;
            G_[idx(i, j)] = v;
            soma += v;
        }

        // Vazamento: y acima do que leva ao fundo da grade.
        // RASTREADO, nunca absorvido no bin de baixo.
        double yl = 1.0 - grid.edge(0)/Ej;
        yl = std::max(0.0, std::min(1.0, yl));
        const double leak = (s_ref > 0.0 && yl < 1.0)
                          ? integra_g(xsec, E_ref, yl, 1.0, abs_piece, quebras)/s_ref
                          : 0.0;
        G_leak_[static_cast<std::size_t>(j)] = leak;
        soma += leak;

        // Identidade discreta: os intervalos particionam [0,1], logo a
        // soma normalizada tem de dar 1.
        residual_[static_cast<std::size_t>(j)] = soma - 1.0;
    }

    const double pior = worst_residual();
    if (!(pior <= consistency_tol)) {
        std::ostringstream m;
        m << "CascadeKernel: identidade de consistencia falhou. "
          << "pior |soma_i G_ij + G_leak_j - 1| = " << pior
          << " > " << consistency_tol
          << ". A cascata nao conservaria numero.";
        throw std::runtime_error(m.str());
    }
}

double CascadeKernel::worst_residual() const
{
    double p = 0.0;
    for (double r : residual_) p = std::max(p, std::fabs(r));
    return p;
}

// =====================================================================
// Z na forma DISCRETA.
//
// Com bins log-uniformes e fluxo em lei de potencia, phi_i (fluxo
// INTEGRADO no bin) satisfaz phi_{i+d}/phi_i = rho^(-d(gamma-1)), e como
// os intervalos de y dependem so de (i-j), a matriz e uma CONVOLUCAO:
// G_ij = G_{j-i}. A lei de potencia e entao AUTOVETOR EXATO do operador
// discreto, com
//
//     sigma_eff,disc = sigma_cc + sigma_nc (1 - Z_disc),
//     Z_disc(gamma)  = soma_d G_d rho^(-d(gamma-1)).
//
// No continuo, (1-y) = E/E' = rho^(-d) no intervalo d, entao Z_disc e
// uma aproximacao de Riemann de
//
//     Z = integral dy g(y) (1-y)^(gamma-1),
//
// avaliando o peso no valor representativo de cada intervalo. Z_disc ->
// Z quando a grade refina; a ordem observada e o que T24 mede.
//
// Isto importa para T22: a solucao da EDO e EXATAMENTE
// exp(-n sigma_eff,disc l), a menos da tolerancia da EDO. Cobrar 1e-10
// contra sigma_eff ANALITICO num numero finito de bins seria cobrar do
// esquema binado algo que ele nao pode dar; contra sigma_eff,disc, sim.
double CascadeKernel::Z_discrete(double gamma) const
{
    // Longe das bordas a convolucao e limpa: tomamos o bin do topo, que
    // e o que tem a cadeia completa de destinos abaixo dele.
    const int j = n_ - 1;
    const double rho = grid_.ratio();

    double z = 0.0;
    for (int i = 0; i <= j; ++i) {
        const double g = G_[idx(i, j)];
        if (g == 0.0) continue;
        const int d = j - i;
        z += g*std::pow(rho, -d*(gamma - 1.0));
    }
    return z;
}



// =====================================================================
namespace {

constexpr double kPiOver2 = 1.5707963267948966;

struct CascadeCtx {
    const CascadeKernel&  ker;
    const Metric&         metric;
    const DensityProfile& profile;
    double r_turn, f_turn, b, sin_i, psi_t, E_inf_ref;
    int    N;
    bool   spherical;
    bool   inbound;      // ramo de entrada: psi = psi_t - psi_off
    double s_hi;
    bool   reversed;     // integra u = s_hi - s

    // grade de E_inf dos bins, so os centros
    std::vector<double> Ec;

    double s_de(double u) const { return reversed ? (s_hi - u) : u; }

    void operator()(double u, const std::vector<double>& y,
                    std::vector<double>& dy) const
    {
        std::fill(dy.begin(), dy.end(), 0.0);

        const double s = s_de(u);
        const double r = r_turn + s*s;

        const RayDerivatives d = ray_derivatives(metric, r, r_turn, b);
        if (!(d.dl_ds > 0.0)) return;

        double theta = kPiOver2;
        if (!spherical) {
            double c = sin_i*std::sin(inbound ? (psi_t - y[static_cast<std::size_t>(N)+1])
                                              : (psi_t + y[static_cast<std::size_t>(N)+1]));
            if (c >  1.0) c =  1.0;
            if (c < -1.0) c = -1.0;
            theta = std::acos(c);
        }

        const double rho = profile.rho(r, theta);
        const double dl  = d.dl_ds;

        // psi sempre avanca (componente N+1)
        dy[static_cast<std::size_t>(N)+1] = d.dpsi_ds;

        if (!(rho > 0.0)) return;

        const double nN    = units::nucleons_per_gram*rho;
        const double boost = 1.0/std::sqrt(metric.f(r));   // E_loc = E_inf * boost

        // coluna (componente N+2)
        dy[static_cast<std::size_t>(N)+2] = rho*dl;

        // perda
        for (int i = 0; i < N; ++i) {
            const double El = Ec[static_cast<std::size_t>(i)]*boost;
            const double st = ker.xsec().sigma_cc(El) + ker.xsec().sigma_nc(El);
            dy[static_cast<std::size_t>(i)] -= nN*st*y[static_cast<std::size_t>(i)]*dl;
        }

        // ganho + vazamento. G_ij nao depende de l: so o coeficiente
        // sigma_nc(E_loc_j) varia.
        for (int j = 0; j < N; ++j) {
            const double pj = y[static_cast<std::size_t>(j)];
            if (pj == 0.0) continue;
            const double El = Ec[static_cast<std::size_t>(j)]*boost;
            const double snc = ker.xsec().sigma_nc(El);
            if (!(snc > 0.0)) continue;

            const double coef = nN*snc*pj*dl;
            for (int i = 0; i <= j; ++i) {
                const double g = ker.G(i, j);
                if (g != 0.0) dy[static_cast<std::size_t>(i)] += coef*g;
            }
            dy[static_cast<std::size_t>(N)] += coef*ker.G_leak(j);
        }
    }
};

} // namespace

CascadeResult transport_cascade(const std::vector<double>& phi0,
                                const CascadeKernel& kernel,
                                const Ray& ray,
                                const Metric& metric,
                                const DensityProfile& profile,
                                const IntegratorOpts& opts)
{
    CascadeResult out;
    const int N = kernel.grid().n();

    if (static_cast<int>(phi0.size()) != N) {
        throw std::runtime_error("transport_cascade: phi0 com tamanho != n_bins");
    }

    out.phi = phi0;

    if (ray.b_cm > 0.0 && metric.is_captured(ray.b_cm)) {
        out.captured = true;
        std::fill(out.phi.begin(), out.phi.end(), 0.0);
        return out;
    }

    const double r_turn = metric.r_turning(ray.b_cm, ray.b_cm);
    const double f_turn = metric.f(r_turn);
    const double r_lo   = std::max(r_turn, profile.r_support_min());
    const double r_hi   = profile.r_support_max();
    if (!(r_hi > r_lo)) return out;    // nao cruza materia

    const double s_lo = std::sqrt(std::max(0.0, r_lo - r_turn));
    const double s_hi = std::sqrt(r_hi - r_turn);

    std::vector<double> Ec(static_cast<std::size_t>(N));
    for (int i = 0; i < N; ++i) Ec[static_cast<std::size_t>(i)] = kernel.grid().center(i);

    // estado: N bins + leakage + psi_off + coluna
    std::vector<double> y(static_cast<std::size_t>(N) + 3, 0.0);
    for (int i = 0; i < N; ++i) y[static_cast<std::size_t>(i)] = phi0[static_cast<std::size_t>(i)];

    // Piso ABSOLUTO por componente, escalado pelo fluxo INICIAL do bin.
    // Ver OdeOpts::abs_tol_per_component para o porque.
    OdeOpts oo;
    oo.rel_tol   = opts.ode_rel_tol;
    oo.max_steps = opts.ode_max_steps;
    oo.abs_tol_per_component.assign(static_cast<std::size_t>(N) + 3, 0.0);
    double phi_max = 0.0;
    for (double v : phi0) phi_max = std::max(phi_max, std::fabs(v));
    for (int i = 0; i < N; ++i) {
        oo.abs_tol_per_component[static_cast<std::size_t>(i)] =
            1.0e-14*std::max(std::fabs(phi0[static_cast<std::size_t>(i)]), 1.0e-30*phi_max);
    }
    oo.abs_tol_per_component[static_cast<std::size_t>(N)]     = 1.0e-14*phi_max;  // leak
    oo.abs_tol_per_component[static_cast<std::size_t>(N)+1]   = 1.0e-12;          // psi
    oo.abs_tol_per_component[static_cast<std::size_t>(N)+2]   = 1.0e-10;          // coluna

    CascadeCtx ctx{kernel, metric, profile,
                   r_turn, f_turn, ray.b_cm, std::sin(ray.inclination_rad),
                   ray.psi_turn_rad, ray.E_inf_GeV, N, profile.is_spherical(),
                   true, s_hi, true, Ec};

    auto roda = [&](double a, double b2) {
        if (!(b2 > a)) return;
        const OdeStats st = dopri54_dyn(
            [&](double u, const std::vector<double>& yy, std::vector<double>& dd){ ctx(u, yy, dd); },
            a, b2, y, oo);
        if (!st.ok) out.tolerance_met = false;
        out.n_evals += st.n_evals;
    };

    // Ramo de ENTRADA: o fluxo chega em s = s_hi e desce ate o ponto de
    // retorno. Integrado em u = s_hi - s, que cresce.
    // psi comeca em psi_t e a EDO acumula psi_off; no ramo de entrada
    // psi = psi_t - psi_off, dai o sinal em CascadeCtx.
    //
    // Aqui psi_off cresce ao descer, ou seja, psi_off(u) e o angulo
    // varrido desde a ENTRADA. Para o mapeamento theta ficar certo,
    // partimos de psi_off = psi_total e descontamos -- por isso o ramo de
    // entrada usa o sinal negativo com psi_off medido do ponto de retorno.
    ctx.inbound = true;  ctx.reversed = true;
    roda(0.0, s_hi - s_lo);

    // No ponto de retorno psi_off volta a zero e o ramo de saida comeca.
    y[static_cast<std::size_t>(N)+1] = 0.0;
    ctx.inbound = false; ctx.reversed = false;
    roda(s_lo, s_hi);

    for (int i = 0; i < N; ++i) out.phi[static_cast<std::size_t>(i)] = y[static_cast<std::size_t>(i)];
    out.leakage_total = y[static_cast<std::size_t>(N)];
    out.column_g_cm2  = y[static_cast<std::size_t>(N)+2];
    return out;
}

} // namespace phasis
