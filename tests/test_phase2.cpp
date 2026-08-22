// =====================================================================
// Fase 2 -- testes de aceitacao T7..T10.
// (T11, nao-regressao, e a suite da Fase 1 rodada sem alteracao.)
// =====================================================================

#include "phasis/cross_section.hpp"
#include "phasis/density.hpp"
#include "phasis/metric.hpp"
#include "phasis/trace.hpp"
#include "phasis/units.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace phasis;

namespace {

constexpr double kPi = 3.141592653589793;

int falhas = 0, total = 0;

void check(const char* nome, double got, double ref, double tol_rel)
{
    ++total;
    const double d = (std::fabs(ref) > 0.0) ? std::fabs(ref) : 1.0;
    const double rel = std::fabs(got - ref)/d;
    const bool ok = (rel <= tol_rel) && std::isfinite(got);
    if (!ok) ++falhas;
    std::printf("    %-44s %-18.11g ref %-18.11g rel %8.2e  %s\n",
                nome, got, ref, rel, ok ? "OK" : "<-- FALHA");
}

void check_bool(const char* nome, bool got, bool ref)
{
    ++total;
    const bool ok = (got == ref);
    if (!ok) ++falhas;
    std::printf("    %-44s %-18s ref %-18s                %s\n",
                nome, got ? "true" : "false", ref ? "true" : "false",
                ok ? "OK" : "<-- FALHA");
}

IntegratorOpts tight(bool defl = false)
{
    IntegratorOpts o;
    o.rel_tol = 1.0e-14;
    o.max_depth = 60;
    o.want_deflection = defl;
    return o;
}

// Ajuste de reta em log-log; devolve a inclinacao.
double slope_loglog(const std::vector<double>& x, const std::vector<double>& y)
{
    const std::size_t n = x.size();
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const double a = std::log(x[i]), b = std::log(y[i]);
        sx += a; sy += b; sxx += a*a; sxy += a*b;
    }
    const double N = static_cast<double>(n);
    return (N*sxy - sx*sy)/(N*sxx - sx*sx);
}

// Schwarzschild forcado a usar a bisecao generica, para valida-la contra
// a forma fechada.
struct SchwarzschildBisect final : Schwarzschild {
    using Schwarzschild::Schwarzschild;
    double r_turning(double b, double rmax) const override {
        return Metric::r_turning(b, rmax);
    }
};

} // namespace

int main()
{
    std::printf("\n=====================================================\n");
    std::printf(" PHASIS -- Fase 2: relatividade geral\n");
    std::printf("=====================================================\n");

    // -----------------------------------------------------------------
    std::printf("\n  T9a  forma fechada de r_t vs bisecao generica\n");
    {
        const double r_s = 1.0e6;
        const Schwarzschild schw(r_s);
        const SchwarzschildBisect bis(r_s);

        for (double frac : {3.0, 10.0, 100.0, 1.0e4}) {
            const double b = frac*r_s;
            char nome[64];
            std::snprintf(nome, sizeof(nome), "r_t  (b/r_s = %g)", frac);
            check(nome, bis.r_turning(b, b), schw.r_turning(b, b), 1.0e-12);
        }
    }

    // -----------------------------------------------------------------
    std::printf("\n  T9b  esfera de fotons e b_crit\n");
    {
        const double r_s = 1.0e6;
        const Schwarzschild schw(r_s);

        const double b_crit_exato = 2.5980762113533160*r_s;
        check("b_crit / r_s", schw.b_crit()/r_s, 2.5980762113533160, 1.0e-15);
        check("r_photon / r_s", schw.r_photon()/r_s, 1.5, 1.0e-15);

        // r_t -> 1.5 r_s quando b -> b_crit+, mas a convergencia e em
        // RAIZ QUADRADA, nao linear: em b = b_crit a cubica tem raiz
        // DUPLA, entao r_t - r_ph ~ sqrt(b - b_crit). Exigir 1e-9 em r_t
        // pediria delta ~ 1e-18 em b, abaixo da dupla precisao. A
        // tolerancia certa e ~sqrt(delta), e o teste forte e a ESCALA.
        std::printf("    aproximacao de b_crit (esperado: raiz dupla, ~sqrt)\n");
        std::printf("    %-14s %-18s %s\n", "delta", "r_t/r_ph - 1", "razao/sqrt(delta)");
        std::vector<double> dd, ee;
        for (int k = 14; k >= 6; k -= 2) {
            const double delta = std::pow(10.0, -k);
            const double rt = schw.r_turning(b_crit_exato*(1.0 + delta), b_crit_exato*2.0);
            const double e = rt/(1.5*r_s) - 1.0;
            std::printf("    %-14.0e %-18.6e %.4f\n", delta, e, e/std::sqrt(delta));
            dd.push_back(delta); ee.push_back(e);
        }
        check("r_t/r_ph - 1 escala como delta^(1/2)", slope_loglog(dd, ee), 0.5, 0.02);
        check("r_t(b_crit*(1+1e-14)) / r_s", 
              schw.r_turning(b_crit_exato*(1.0 + 1.0e-14), b_crit_exato*2.0)/r_s,
              1.5, 1.0e-6);

        // transicao exata, 10 digitos
        check_bool("capturado em b_crit*(1 - 1e-12)",
                   schw.is_captured(b_crit_exato*(1.0 - 1.0e-12)), true);
        check_bool("livre     em b_crit*(1 + 1e-12)",
                   schw.is_captured(b_crit_exato*(1.0 + 1.0e-12)), false);

        // bisecao por busca: mesma transicao
        double lo = 2.0*r_s, hi = 4.0*r_s;
        for (int i = 0; i < 200; ++i) {
            const double m = 0.5*(lo + hi);
            if (schw.is_captured(m)) lo = m; else hi = m;
        }
        check("transicao localizada / r_s", 0.5*(lo + hi)/r_s,
              2.5980762113533160, 1.0e-10);

        // T(r_t) = 2 r_t - 3 r_s, e zera na esfera de fotons
        check("T(r_t) em b_crit  (deve ser ~0)",
              schw.turning_factor_at_turn(1.5*r_s)/r_s, 0.0, 1.0e-15);
        check("T(r_t) em r_t = 5 r_s",
              schw.turning_factor_at_turn(5.0*r_s)/r_s, 7.0, 1.0e-15);

        // Result: captura devolve tau infinito, P = 0, sem NaN
        const UniformBall ball(1.0, 50.0*r_s);
        const PowerLawCrossSection xs(1.0e-33, 1.0, 0.0);
        Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = 2.0*r_s;
        const Result r = trace_ray(ray, schw, ball, xs, tight());
        check_bool("Result.captured", r.captured, true);
        check_bool("tau = +inf", std::isinf(r.tau) && r.tau > 0.0, true);
        check("P_surv", r.P_surv, 0.0, 0.0);
    }

    // -----------------------------------------------------------------
    std::printf("\n  T9c  regime quase-critico: divergencia logaritmica\n");
    {
        const double r_s = 1.0e6;
        const Schwarzschild schw(r_s);
        const double b_c = schw.b_crit();

        const UniformBall ball(1.0e-3, 1.0e4*r_s);
        const PowerLawCrossSection xs(1.0e-30, 1.0, 0.0);

        // O observavel limpo da divergencia e o ENROLAMENTO: geometria
        // pura, sem depender do perfil de materia. Delta_phi diverge
        // logaritmicamente em delta = b/b_crit - 1, porque o raio da
        // cada vez mais voltas na esfera de fotons antes de escapar.
        std::printf("    %-12s %-14s %-12s %-10s %-9s %s\n",
                    "b/b_c - 1", "tau", "T(r_t)/r_t", "voltas", "near_crit", "finito");

        std::vector<double> lnd, voltas;
        for (int k = 4; k <= 14; k += 2) {
            const double delta = std::pow(10.0, -k);
            Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = b_c*(1.0 + delta);

            // Tolerancia folgada de proposito: no regime quase-critico o
            // integrando tem estrutura em varias escalas, e pedir 1e-14
            // e pedir abaixo do epsilon de maquina.
            IntegratorOpts o;
            o.rel_tol = 1.0e-10;
            o.max_depth = 40;
            o.want_deflection = true;

            const Result r = trace_ray(ray, schw, ball, xs, o);
            const double eps = schw.turning_factor_at_turn(r.r_min_cm)/r.r_min_cm;

            const bool fin = std::isfinite(r.tau) && std::isfinite(r.deflection_rad)
                          && r.tau > 0.0;
            std::printf("    %-12.0e %-14.6f %-12.3e %-10.4f %-9s %s\n",
                        delta, r.tau, eps, r.winding_turns,
                        r.near_critical ? "sim" : "nao", fin ? "sim" : "NAO");

            ++total; if (!fin) { ++falhas; std::printf("      <-- FALHA\n"); }

            lnd.push_back(-std::log(delta));
            voltas.push_back(r.winding_turns);
        }

        // Delta_phi = a + b*(-ln delta). Ajuste linear + R^2.
        const std::size_t n = lnd.size();
        double sx=0, sy=0, sxx=0, sxy=0;
        for (std::size_t i = 0; i < n; ++i) {
            sx += lnd[i]; sy += voltas[i]; sxx += lnd[i]*lnd[i]; sxy += lnd[i]*voltas[i];
        }
        const double N = static_cast<double>(n);
        const double incl = (N*sxy - sx*sy)/(N*sxx - sx*sx);
        const double mx = sx/N, my = sy/N;
        double ss_tot=0, ss_res=0;
        for (std::size_t i = 0; i < n; ++i) {
            const double pred = my + incl*(lnd[i] - mx);
            ss_res += (voltas[i]-pred)*(voltas[i]-pred);
            ss_tot += (voltas[i]-my)*(voltas[i]-my);
        }
        std::printf("\n    voltas = a + b*(-ln delta):  b = %.6f,  R^2 = %.9f\n",
                    incl, 1.0 - ss_res/ss_tot);

        // Coeficiente analitico, e um teste forte de verdade.
        //
        // No limite de deflexao forte (Bozza), Delta_phi = -a ln(b/b_c - 1) + c,
        // e para Schwarzschild o coeficiente vale exatamente
        //
        //     a = 1
        //
        // (e o expoente de Lyapunov da orbita circular instavel na
        // esfera de fotons). Logo voltas = Delta_phi/(2pi) tem
        // inclinacao 1/(2 pi) contra -ln(delta).
        check("linearidade das voltas em -ln(delta)  (R^2)",
              1.0 - ss_res/ss_tot, 1.0, 1.0e-6);
        check("inclinacao das voltas = 1/(2 pi)  [Bozza: a = 1]",
              incl, 1.0/(2.0*kPi), 1.0e-3);

        // captura estrita logo abaixo do critico
        Ray rc; rc.E_inf_GeV = 1.0e9; rc.b_cm = b_c*(1.0 - 1.0e-10);
        const Result cap = trace_ray(rc, schw, ball, xs, tight());
        check_bool("b < b_crit: captured", cap.captured, true);
        check_bool("b < b_crit: tau = +inf, P = 0, sem NaN",
                   std::isinf(cap.tau) && cap.P_surv == 0.0 && !std::isnan(cap.P_surv),
                   true);
    }

    // -----------------------------------------------------------------
    std::printf("\n  T7  convergencia newtoniana: desvio ~ (r_s/r_t)^1\n");
    {
        const double R = 1.0e12, rho0 = 1.0, sig = 1.0e-33;
        const double b = 0.5*R;
        const UniformBall ball(rho0, R);
        const PowerLawCrossSection xs(sig, 1.0, 0.0);   // sigma constante
        const Minkowski flat;

        Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = b;
        const double tau_M = trace_ray(ray, flat, ball, xs, tight()).tau;

        std::vector<double> xs_fit, ys_fit;
        std::printf("    %-14s %-16s %s\n", "r_s/r_t", "desvio relativo", "");
        for (int k = 12; k >= 3; --k) {
            const double ratio = std::pow(10.0, -k);
            const Schwarzschild schw(ratio*b);
            const double tau_S = trace_ray(ray, schw, ball, xs, tight()).tau;
            const double dev = std::fabs(tau_S - tau_M)/tau_M;
            std::printf("    %-14.1e %-16.6e\n", ratio, dev);
            if (k <= 9) { xs_fit.push_back(ratio); ys_fit.push_back(dev); }
        }
        const double m = slope_loglog(xs_fit, ys_fit);
        std::printf("\n");
        check("inclinacao log-log (decadas 1e-9..1e-3)", m, 1.0, 0.02);
    }

    // -----------------------------------------------------------------
    std::printf("\n  T8  deflexao de campo fraco: Delta_phi -> 2 r_s / b\n");
    {
        const double r_s = 1.0e5;
        const Schwarzschild schw(r_s);
        const UniformBall vazio(1.0, 1.0e-3*r_s);   // nenhum cruzamento
        const PowerLawCrossSection xs(1.0e-33, 1.0, 0.0);

        std::vector<double> xs_fit, ys_fit;
        std::printf("    %-12s %-18s %-18s %s\n",
                    "b/r_s", "Delta_phi", "2 r_s / b", "erro rel");
        for (double frac : {1.0e3, 1.0e4, 1.0e5}) {
            const double b = frac*r_s;
            Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = b;
            const Result r = trace_ray(ray, schw, vazio, xs, tight(true));

            const double esperado = 2.0*r_s/b;
            const double err = std::fabs(r.deflection_rad - esperado)/esperado;
            std::printf("    %-12.0e %-18.10e %-18.10e %.3e\n",
                        frac, r.deflection_rad, esperado, err);

            char nome[64];
            std::snprintf(nome, sizeof(nome), "Delta_phi  (b/r_s = %.0e)", frac);
            check(nome, r.deflection_rad, esperado, 0.01);

            xs_fit.push_back(1.0/frac);
            ys_fit.push_back(err);
        }
        const double m = slope_loglog(xs_fit, ys_fit);
        std::printf("\n");
        check("erro cai como (r_s/b)^1", m, 1.0, 0.05);

        // Teste mais forte que o pedido: nao so o termo dominante, mas o
        // COEFICIENTE DE SEGUNDA ORDEM. A expansao no parametro de
        // impacto (nao em r_min) e
        //
        //     Delta_phi = 2 (r_s/b) + (15 pi/16) (r_s/b)^2 + ...
        //
        // Extraindo o coeficiente do residuo: se ele bate, a geometria
        // esta certa a segunda ordem, nao so no limite.
        std::printf("    coeficiente de 2a ordem  (esperado 15 pi/16 = %.6f)\n",
                    15.0*kPi/16.0);
        for (double frac : {1.0e3, 1.0e4, 1.0e5}) {
            const double b = frac*r_s;
            Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = b;
            const double d = trace_ray(ray, schw, vazio, xs, tight(true)).deflection_rad;
            const double x = r_s/b;
            const double c2 = (d - 2.0*x)/(x*x);
            std::printf("    b/r_s = %-10.0e c2 = %.6f\n", frac, c2);
            if (frac >= 1.0e4) {
                char nome[64];
                std::snprintf(nome, sizeof(nome), "c2  (b/r_s = %.0e)", frac);
                check(nome, c2, 15.0*kPi/16.0, 2.0e-3);
            }
        }

        // Minkowski: deflexao identicamente zero
        const Minkowski flat;
        Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = 1.0e8;
        const Result rm = trace_ray(ray, flat, vazio, xs, tight(true));
        check("deflexao em Minkowski", rm.deflection_rad, 0.0, 1.0e-12);
    }

    // -----------------------------------------------------------------
    std::printf("\n  T10  redshift: sigma DENTRO da integral\n");
    {
        const double R = 1.0e12, rho0 = 1.0;
        const double sigma0 = 1.0e-33, E0 = 1.0e6, alpha = 0.4;
        const double E_inf = 1.0e9;
        const double b = 0.5*R;

        const UniformBall ball(rho0, R);
        const PowerLawCrossSection xs_E(sigma0, E0, alpha);

        // O jeito ERRADO: sigma congelada em E_inf, fatorada para fora.
        const PowerLawCrossSection xs_fixo(xs_E.sigma_tot(E_inf), 1.0, 0.0);

        Ray ray; ray.E_inf_GeV = E_inf; ray.b_cm = b;

        std::vector<double> xs_fit, ys_fit;
        std::printf("    %-14s %-16s %s\n", "r_s/r_t", "tau_certo/tau_errado - 1", "");
        for (int k = 8; k >= 3; --k) {
            const double ratio = std::pow(10.0, -k);
            const Schwarzschild schw(ratio*b);

            const double tau_certo  = trace_ray(ray, schw, ball, xs_E,   tight()).tau;
            const double tau_errado = trace_ray(ray, schw, ball, xs_fixo, tight()).tau;
            const double d = tau_certo/tau_errado - 1.0;
            std::printf("    %-14.1e %-16.6e\n", ratio, d);
            xs_fit.push_back(ratio);
            ys_fit.push_back(std::fabs(d));
        }
        std::printf("\n");

        // (i) a diferenca existe
        check_bool("diferenca nao nula em r_s/r_t = 1e-3", ys_fit.back() > 1.0e-6, true);
        // (ii) e cresce linearmente com r_s/r_t
        check("inclinacao log-log", slope_loglog(xs_fit, ys_fit), 1.0, 0.02);
        // (iii) sinal: E_loc > E_inf e sigma cresce com E, logo tau_certo > tau_errado
        const Schwarzschild schw(1.0e-3*b);
        const double tc = trace_ray(ray, schw, ball, xs_E,   tight()).tau;
        const double te = trace_ray(ray, schw, ball, xs_fixo, tight()).tau;
        check_bool("tau_certo > tau_errado", tc > te, true);

        // (iv) coeficiente: tau_certo/tau_errado - 1 ~ (alpha/2) <r_s/r>
        //      com r entre r_t e R, o valor medio fica entre alpha/2 * r_s/R
        //      e alpha/2 * r_s/r_t. Checagem de ordem de grandeza.
        const double d = tc/te - 1.0;
        const double lo = 0.5*alpha*(1.0e-3*b)/R;
        const double hi = 0.5*alpha*(1.0e-3*b)/b;
        check_bool("coeficiente entre (alpha/2)r_s/R e (alpha/2)r_s/r_t",
                   d > lo && d < hi, true);
    }

    std::printf("\n=====================================================\n");
    std::printf("  %d verificacoes, %d falhas\n", total, falhas);
    std::printf("  RESULTADO: %s\n", falhas == 0 ? "OK" : "FALHA");
    std::printf("=====================================================\n\n");
    return falhas == 0 ? 0 : 1;
}
