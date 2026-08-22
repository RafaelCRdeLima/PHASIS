#include "phasis/cascade.hpp"
#include "phasis/integrate.hpp"
#include "phasis/ode.hpp"
#include "phasis/trace.hpp"
#include "phasis/units.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
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
                             int /*n_y_quad, obsoleto: a quadratura e adaptativa*/,
                             int n_E_nodes,
                             double E_lo, double E_hi)
    : grid_(grid), xsec_(xsec), n_(grid.n())
{
    // ---- nos em E ----------------------------------------------------
    //
    // Forma independente de E: um no so, e G(i,j) vale em toda parte.
    // Forma dependente de E: nos log-espacados, e G(i,j,E) interpola
    // LINEARMENTE em ln E -- ver o comentario da classe sobre por que a
    // interpolacao tem de ser linear.
    if (xsec.shape_is_E_independent()) {
        nE_ = 1;
        E_nodes_.assign(1, grid.center(0));
    } else {
        if (n_E_nodes < 2) n_E_nodes = 2;
        // Faixa automatica: a grade toda, com folga de 2x para cima e
        // para baixo. O boost maximo do redshift em Schwarzschild e
        // sqrt(3) = 1.73, entao 2 cobre com margem.
        double lo = (E_lo > 0.0) ? E_lo : 0.5*grid.E_min();
        double hi = (E_hi > 0.0) ? E_hi : 2.0*grid.E_max();
        lo = std::max(lo, xsec.E_domain_min());
        hi = std::min(hi, xsec.E_domain_max());
        if (!(hi > lo)) {
            std::ostringstream m;
            m << "CascadeKernel: a faixa de E pedida [" << lo << ", " << hi
              << "] esta vazia depois de cortada pelo dominio de "
              << xsec.name() << ". A tabela nao cobre a grade de energia.";
            throw std::runtime_error(m.str());
        }
        nE_ = n_E_nodes;
        E_nodes_.resize(static_cast<std::size_t>(nE_));
        for (int e = 0; e < nE_; ++e) {
            E_nodes_[static_cast<std::size_t>(e)] =
                lo*std::pow(hi/lo, static_cast<double>(e)/(nE_ - 1));
        }
    }
    lnE_.resize(E_nodes_.size());
    for (std::size_t e = 0; e < E_nodes_.size(); ++e) lnE_[e] = std::log(E_nodes_[e]);

    const std::size_t NN = static_cast<std::size_t>(n_)*static_cast<std::size_t>(n_);
    G_.assign(static_cast<std::size_t>(nE_)*NN, 0.0);
    G_leak_.assign(static_cast<std::size_t>(nE_)*static_cast<std::size_t>(n_), 0.0);
    residual_.assign(static_cast<std::size_t>(nE_)*static_cast<std::size_t>(n_), 0.0);

    // Caso degenerado: sem NC nao ha o que redistribuir. A identidade
    // NORMALIZADA (soma G = 1) e 0/0 e nao faz sentido; a identidade
    // real, soma_i K_ij + leak_j = sigma_nc = 0, vale trivialmente.
    // G fica tudo zero e o residuo e zero, nao -1.
    if (!(xsec.sigma_nc(E_nodes_[0]) > 0.0) && nE_ == 1) {
        no_regeneration_ = true;
        return;
    }

    const std::vector<double> quebras = xsec.y_breakpoints();

    for (int e = 0; e < nE_; ++e) {
        const double E_ref = E_nodes_[static_cast<std::size_t>(e)];
        const double s_ref = xsec.sigma_nc(E_ref);
        if (!(s_ref > 0.0)) {
            if (nE_ == 1) { no_regeneration_ = true; return; }
            throw std::runtime_error(
                "CascadeKernel: sigma_nc <= 0 num no interno da tabela de E.");
        }

        // Tolerancia ABSOLUTA por pedaco, escalada por sigma_nc. E o que
        // garante que a soma de muitos pedacos ainda feche a identidade.
        const double abs_piece = 1.0e-16*std::fabs(s_ref);

        for (int j = 0; j < n_; ++j) {
            const double Ej_inf = grid.center(j);

            // Os intervalos em y NAO dependem de E_loc: y e invariante
            // sob redshift e as bordas dos bins estao em E_inf.
            double soma = 0.0;
            for (int i = 0; i <= j; ++i) {
                double y0 = 1.0 - grid.edge(i + 1)/Ej_inf;
                double y1 = 1.0 - grid.edge(i)/Ej_inf;
                y0 = std::max(0.0, y0);
                y1 = std::min(1.0, y1);
                if (!(y1 > y0)) continue;

                const double v = integra_g(xsec, E_ref, y0, y1, abs_piece, quebras);
                G_[idx(e, i, j)] = v;
                soma += v;
            }

            // Vazamento: y acima do que leva ao fundo da grade.
            // RASTREADO, nunca absorvido no bin de baixo.
            double yl = 1.0 - grid.edge(0)/Ej_inf;
            yl = std::max(0.0, std::min(1.0, yl));
            const double leak = (yl < 1.0)
                ? integra_g(xsec, E_ref, yl, 1.0, abs_piece, quebras)
                : 0.0;
            G_leak_[static_cast<std::size_t>(e)*static_cast<std::size_t>(n_)
                    + static_cast<std::size_t>(j)] = leak;
            soma += leak;

            // NORMALIZACAO PELA SOMA DOS PEDACOS, nao por sigma_nc.
            //
            // Os intervalos particionam [0,1] exatamente, entao a soma
            // dos pedacos E a integral de dsigma/dy em [0,1] -- calculada
            // com o MESMO integrador adaptativo que calculou cada pedaco.
            // Dividir por sigma_nc, que num kernel de tabela vem de um
            // trapezio sobre os nos em y, misturaria dois integradores e
            // a identidade nao fecharia: sobraria a diferenca entre eles,
            // que aparece como fluxo criado ou destruido ao longo do raio.
            //
            // A diferenca entre os dois NAO e escondida: fica em
            // norm_mismatch(), e mede a resolucao da grade em y.
            if (!(soma > 0.0)) {
                residual_[static_cast<std::size_t>(e)*static_cast<std::size_t>(n_)
                          + static_cast<std::size_t>(j)] = 0.0;
                continue;
            }
            for (int i = 0; i <= j; ++i) G_[idx(e, i, j)] /= soma;
            G_leak_[static_cast<std::size_t>(e)*static_cast<std::size_t>(n_)
                    + static_cast<std::size_t>(j)] /= soma;

            double conf = 0.0;
            for (int i = 0; i <= j; ++i) conf += G_[idx(e, i, j)];
            conf += G_leak_[static_cast<std::size_t>(e)*static_cast<std::size_t>(n_)
                            + static_cast<std::size_t>(j)];
            residual_[static_cast<std::size_t>(e)*static_cast<std::size_t>(n_)
                      + static_cast<std::size_t>(j)] = conf - 1.0;

            if (j == n_ - 1) {
                norm_mismatch_.push_back(soma/s_ref - 1.0);
            }
        }
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

void CascadeKernel::locate(double E, int& e0, int& e1, double& t) const
{
    if (nE_ == 1) { e0 = e1 = 0; t = 0.0; return; }

    const double x = std::log(E);
    if (!(x >= lnE_.front()) || !(x <= lnE_.back())) {
        std::ostringstream m;
        m << "CascadeKernel: E = " << E << " GeV fora da faixa tabelada ["
          << std::exp(lnE_.front()) << ", " << std::exp(lnE_.back())
          << "]. O kernel nao extrapola: a forma em y fora da faixa nao foi "
             "medida, e inventa-la produziria uma cascata plausivel e errada.";
        throw std::out_of_range(m.str());
    }
    const auto it = std::upper_bound(lnE_.begin(), lnE_.end(), x);
    std::size_t k = static_cast<std::size_t>(it - lnE_.begin());
    if (k == 0) k = 1;
    if (k >= lnE_.size()) k = lnE_.size() - 1;
    e0 = static_cast<int>(k) - 1;
    e1 = static_cast<int>(k);
    t  = (x - lnE_[static_cast<std::size_t>(e0)])
       / (lnE_[static_cast<std::size_t>(e1)] - lnE_[static_cast<std::size_t>(e0)]);
}

double CascadeKernel::G(int i, int j) const
{
    if (nE_ > 1) {
        throw std::logic_error(
            "CascadeKernel::G(i,j) sem energia: a forma em y desta secao de "
            "choque DEPENDE de E. Use G(i,j,E). Devolver o primeiro no seria "
            "dar um numero plausivel e errado.");
    }
    return G_[idx(0, i, j)];
}

double CascadeKernel::G_leak(int j) const
{
    if (nE_ > 1) {
        throw std::logic_error(
            "CascadeKernel::G_leak(j) sem energia: a forma depende de E. "
            "Use G_leak(j,E).");
    }
    return G_leak_[static_cast<std::size_t>(j)];
}

double CascadeKernel::G(int i, int j, double E) const
{
    int e0, e1; double t;
    locate(E, e0, e1, t);
    const double a = G_[idx(e0, i, j)];
    if (e0 == e1) return a;
    return a + t*(G_[idx(e1, i, j)] - a);
}

double CascadeKernel::G_leak(int j, double E) const
{
    int e0, e1; double t;
    locate(E, e0, e1, t);
    const std::size_t nn = static_cast<std::size_t>(n_);
    const double a = G_leak_[static_cast<std::size_t>(e0)*nn + static_cast<std::size_t>(j)];
    if (e0 == e1) return a;
    const double b = G_leak_[static_cast<std::size_t>(e1)*nn + static_cast<std::size_t>(j)];
    return a + t*(b - a);
}

double CascadeKernel::worst_norm_mismatch() const
{
    double p = 0.0;
    for (double v : norm_mismatch_) p = std::max(p, std::fabs(v));
    return p;
}

double CascadeKernel::consistency_residual(int j) const
{
    double p = 0.0;
    const std::size_t nn = static_cast<std::size_t>(n_);
    for (int e = 0; e < nE_; ++e) {
        const double r = residual_[static_cast<std::size_t>(e)*nn
                                   + static_cast<std::size_t>(j)];
        if (std::fabs(r) > std::fabs(p)) p = r;
    }
    return p;
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
    if (nE_ > 1) {
        throw std::logic_error(
            "CascadeKernel::Z_discrete: so faz sentido com forma independente "
            "de E; com forma dependente de E nao ha um unico Z.");
    }
    // Longe das bordas a convolucao e limpa: tomamos o bin do topo, que
    // e o que tem a cadeia completa de destinos abaixo dele.
    const int j = n_ - 1;
    const double rho = grid_.ratio();

    double z = 0.0;
    for (int i = 0; i <= j; ++i) {
        const double g = G_[idx(0, i, j)];
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
    bool   freeze;       // T25: usa E_inf no lugar de E_loc

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
        const double boost = freeze ? 1.0 : 1.0/std::sqrt(metric.f(r));

        // coluna (componente N+2)
        dy[static_cast<std::size_t>(N)+2] = rho*dl;

        // perda
        for (int i = 0; i < N; ++i) {
            const double El = Ec[static_cast<std::size_t>(i)]*boost;
            const double st = ker.xsec().sigma_cc(El) + ker.xsec().sigma_nc(El);
            dy[static_cast<std::size_t>(i)] -= nN*st*y[static_cast<std::size_t>(i)]*dl;
        }

        // Ganho + vazamento.
        //
        // Os INTERVALOS em y que ligam j a i nao dependem de l -- y e
        // invariante sob redshift e as bordas dos bins estao em E_inf.
        // Com kernel analitico a FORMA tambem nao depende, e G_ij e uma
        // constante. Com tabela real a forma depende de E, e G_ij(E_loc)
        // e interpolado em ln E a partir dos nos montados uma vez na
        // construcao -- nunca reintegrado por passo.
        for (int j = 0; j < N; ++j) {
            const double pj = y[static_cast<std::size_t>(j)];
            if (pj == 0.0) continue;
            const double El = Ec[static_cast<std::size_t>(j)]*boost;
            const double snc = ker.xsec().sigma_nc(El);
            if (!(snc > 0.0)) continue;

            const double coef = nN*snc*pj*dl;
            for (int i = 0; i <= j; ++i) {
                const double g = ker.G(i, j, El);
                if (g != 0.0) dy[static_cast<std::size_t>(i)] += coef*g;
            }
            dy[static_cast<std::size_t>(N)] += coef*ker.G_leak(j, El);
        }
    }
};

} // namespace

CascadeResult transport_cascade(const std::vector<double>& phi0,
                                const CascadeKernel& kernel,
                                const Ray& ray,
                                const Metric& metric,
                                const DensityProfile& profile,
                                const IntegratorOpts& opts,
                                bool freeze_redshift)
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
                   true, s_hi, true, freeze_redshift, Ec};

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


// =====================================================================
double z_discrete_row(const EnergyGrid& grid,
                      const DifferentialCrossSection& xsec,
                      double gamma)
{
    const int j = grid.n() - 1;
    const double Ej   = grid.center(j);
    const double rho  = grid.ratio();
    const double Eref = grid.center(0);
    const double sref = xsec.sigma_nc(Eref);
    if (!(sref > 0.0)) return 0.0;

    const std::vector<double> quebras = xsec.y_breakpoints();
    const double abs_piece = 1.0e-16*std::fabs(sref);

    double z = 0.0;
    for (int i = 0; i <= j; ++i) {
        double y0 = std::max(0.0, 1.0 - grid.edge(i + 1)/Ej);
        double y1 = std::min(1.0, 1.0 - grid.edge(i)/Ej);
        if (!(y1 > y0)) continue;
        const double g = integra_g(xsec, Eref, y0, y1, abs_piece, quebras)/sref;
        if (g == 0.0) continue;
        z += g*std::pow(rho, -(j - i)*(gamma - 1.0));
    }
    return z;
}


// =====================================================================
namespace {

const char* kChavesObrigatorias[] = {
    "convention_y", "target", "projectile", "current",
    "units_sigma", "units_E", "M_Z_GeV", "dipole_model", "generated_by"
};

} // namespace

void write_dsigma_table(const std::string& path,
                        const DifferentialCrossSection& xsec,
                        double E_min, double E_max, int nE,
                        double y_min, int nY,
                        const std::string& gerado_por)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("write_dsigma_table: nao consegui escrever " + path);

    out << "# convention_y = (E_in - E_out)/E_in\n"
        << "# target       = isoscalar_nucleon\n"
        << "# projectile   = nu\n"
        << "# current      = NC\n"
        << "# units_sigma  = cm^2\n"
        << "# units_E      = GeV\n"
        << "# M_Z_GeV      = 91.1876\n"
        << "# dipole_model = " << xsec.name() << "\n"
        << "# generated_by = " << gerado_por << "\n"
        << "# E_GeV y dsigma_dy_cm2\n";
    out.precision(17);

    for (int a = 0; a < nE; ++a) {
        const double E = E_min*std::pow(E_max/E_min, static_cast<double>(a)/(nE-1));
        for (int b = 0; b < nY; ++b) {
            const double y = y_min*std::pow(1.0/y_min, static_cast<double>(b)/(nY-1));
            out << E << " " << y << " " << xsec.dsigma_nc_dy(E, y) << "\n";
        }
    }
}

TableDifferentialCrossSection::TableDifferentialCrossSection(
    const std::string& path, const Expect& esperado,
    std::shared_ptr<const CrossSection> sigma_cc)
    : path_(path), cc_(std::move(sigma_cc))
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("TableDifferential: nao consegui abrir " + path);

    std::vector<double> E, y, d;
    std::string linha;
    while (std::getline(in, linha)) {
        if (!linha.empty() && linha[0] == '#') {
            const auto eq = linha.find('=');
            if (eq == std::string::npos) continue;
            auto trim = [](std::string t){
                const auto a = t.find_first_not_of(" \t#\r\n");
                if (a == std::string::npos) return std::string{};
                const auto b = t.find_last_not_of(" \t\r\n");
                return t.substr(a, b - a + 1); };
            meta_.emplace_back(trim(linha.substr(0, eq)), trim(linha.substr(eq + 1)));
            continue;
        }
        std::istringstream ss(linha);
        double a, b, c;
        if (ss >> a >> b >> c) { E.push_back(a); y.push_back(b); d.push_back(c); }
    }

    // Metadados: ausencia de QUALQUER chave e excecao, nunca default.
    for (const char* ch : kChavesObrigatorias) {
        bool achou = false;
        for (const auto& kv : meta_) if (kv.first == ch) achou = true;
        if (!achou) {
            throw std::runtime_error(
                std::string("TableDifferential: falta o metadado obrigatorio '") + ch
                + "' em " + path + ". O formato exige a convencao explicita: "
                "e ai que mora o fator 2.");
        }
    }
    if (meta("convention_y") != esperado.convention_y) {
        throw std::runtime_error("TableDifferential: convention_y = '" + meta("convention_y")
            + "' mas o chamador espera '" + esperado.convention_y + "'");
    }
    if (meta("current") != esperado.current) {
        throw std::runtime_error("TableDifferential: current = '" + meta("current")
            + "' mas o chamador espera '" + esperado.current + "'");
    }

    // grade: E varia devagar, y rapido (ordem de escrita)
    std::vector<double> Eu, yu;
    for (double v : E) if (Eu.empty() || v != Eu.back()) Eu.push_back(v);
    for (double v : y) { if (yu.size() && v == yu.front()) break; yu.push_back(v); }
    const std::size_t nE = Eu.size(), nY = yu.size();
    if (nE < 2 || nY < 2 || d.size() != nE*nY) {
        throw std::runtime_error("TableDifferential: grade inconsistente em " + path);
    }

    lnE_.reserve(nE); for (double v : Eu) lnE_.push_back(std::log(v));
    lny_.reserve(nY); for (double v : yu) lny_.push_back(std::log(v));
    lnD_.assign(nE*nY, -1.0e300);
    for (std::size_t q = 0; q < d.size(); ++q)
        if (d[q] > 0.0) lnD_[q] = std::log(d[q]);

    // sigma_nc por no de E: integral em y da propria tabela, em log y
    snc_.assign(nE, 0.0);
    for (std::size_t a = 0; a < nE; ++a) {
        double acc = 0.0;
        for (std::size_t b = 0; b + 1 < nY; ++b) {
            const double d0 = d[a*nY + b], d1 = d[a*nY + b + 1];
            acc += 0.5*(d0 + d1)*(yu[b+1] - yu[b]);
        }
        snc_[a] = acc;
    }
}

std::vector<double> TableDifferentialCrossSection::y_breakpoints() const
{
    // So a BORDA DO SUPORTE. Abaixo de y_min a tabela devolve zero e
    // acima dela o valor e finito: isso e um SALTO, e um salto dentro de
    // um painel de Simpson nao cancela entre intervalos vizinhos -- e
    // exatamente a classe de bug que custou um fator 670 de ruido no
    // dipole.
    //
    // Os nos internos da grade em y sao apenas quebras de DERIVADA (a
    // interpolacao bilinear e C^0). Simpson sobre uma quebra de derivada
    // erra em O(h^3), nao em O(h): devolve-los todos aqui subdividiria
    // cada um dos ~1800 intervalos em 120 pedacos, por um erro que ja e
    // menor que o do proprio trapezio que gerou sigma_nc.
    return { std::exp(lny_.front()), std::exp(lny_.back()) };
}

const std::string& TableDifferentialCrossSection::meta(const std::string& chave) const
{
    for (const auto& kv : meta_) if (kv.first == chave) return kv.second;
    throw std::runtime_error("TableDifferential: metadado ausente: " + chave);
}

double TableDifferentialCrossSection::dsigma_nc_dy(double E_GeV, double y) const
{
    if (!(E_GeV > 0.0) || !(y > 0.0)) return 0.0;
    const double x = std::log(E_GeV), t = std::log(y);

    // Sem extrapolacao, pelo mesmo motivo da tabela de CC.
    if (x < lnE_.front() || x > lnE_.back() ||
        t < lny_.front() || t > lny_.back()) {
        return 0.0;
    }

    const std::size_t nY = lny_.size();
    auto ia = static_cast<std::size_t>(
        std::upper_bound(lnE_.begin(), lnE_.end(), x) - lnE_.begin());
    if (ia == 0) ia = 1;
    if (ia >= lnE_.size()) ia = lnE_.size() - 1;
    auto ib = static_cast<std::size_t>(
        std::upper_bound(lny_.begin(), lny_.end(), t) - lny_.begin());
    if (ib == 0) ib = 1;
    if (ib >= nY) ib = nY - 1;

    const std::size_t a0 = ia - 1, b0 = ib - 1;
    const double u = (x - lnE_[a0])/(lnE_[ia] - lnE_[a0]);
    const double v = (t - lny_[b0])/(lny_[ib] - lny_[b0]);

    const double f00 = lnD_[a0*nY + b0], f01 = lnD_[a0*nY + ib];
    const double f10 = lnD_[ia*nY + b0], f11 = lnD_[ia*nY + ib];
    if (f00 < -1.0e299 || f01 < -1.0e299 || f10 < -1.0e299 || f11 < -1.0e299) return 0.0;

    return std::exp((1-u)*(1-v)*f00 + (1-u)*v*f01 + u*(1-v)*f10 + u*v*f11);
}

double TableDifferentialCrossSection::sigma_nc(double E_GeV) const
{
    if (!(E_GeV > 0.0)) return 0.0;
    const double x = std::log(E_GeV);
    if (x < lnE_.front() || x > lnE_.back()) return 0.0;
    auto ia = static_cast<std::size_t>(
        std::upper_bound(lnE_.begin(), lnE_.end(), x) - lnE_.begin());
    if (ia == 0) ia = 1;
    if (ia >= lnE_.size()) ia = lnE_.size() - 1;
    const std::size_t a0 = ia - 1;
    const double u = (x - lnE_[a0])/(lnE_[ia] - lnE_[a0]);
    return std::exp((1-u)*std::log(snc_[a0]) + u*std::log(snc_[ia]));
}

double TableDifferentialCrossSection::sigma_cc(double E_GeV) const
{
    return cc_ ? cc_->sigma_tot(E_GeV) : 0.0;
}

} // namespace phasis
