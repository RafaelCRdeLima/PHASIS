// =====================================================================
// Fase 1 -- testes de aceitacao T1..T6.
//
// Framework minimo proprio, sem dependencia externa.
// =====================================================================

#include "phasis/cross_section.hpp"
#include "phasis/density.hpp"
#include "phasis/metric.hpp"
#include "phasis/trace.hpp"
#include "phasis/units.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

using namespace phasis;

// ---------------------------------------------------------------------
namespace {

int falhas = 0;
int total  = 0;

void check(const std::string& nome, double got, double esperado, double tol_rel)
{
    ++total;
    const double denom = (std::fabs(esperado) > 0.0) ? std::fabs(esperado) : 1.0;
    const double rel = std::fabs(got - esperado)/denom;
    const bool ok = (rel <= tol_rel) && std::isfinite(got);
    if (!ok) ++falhas;
    std::printf("    %-46s %-16.10g  ref %-16.10g  rel %8.2e  %s\n",
                nome.c_str(), got, esperado, rel, ok ? "OK" : "<-- FALHA");
}

void check_bool(const std::string& nome, bool got, bool esperado)
{
    ++total;
    const bool ok = (got == esperado);
    if (!ok) ++falhas;
    std::printf("    %-46s %-16s  ref %-16s              %s\n",
                nome.c_str(), got ? "true" : "false", esperado ? "true" : "false",
                ok ? "OK" : "<-- FALHA");
}

IntegratorOpts opts_tight()
{
    IntegratorOpts o;
    o.rel_tol   = 1.0e-13;
    o.max_depth = 60;
    return o;
}

} // namespace

// ---------------------------------------------------------------------
int main()
{
    const Minkowski flat;
    const double N_A = units::N_A;

    std::printf("\n=====================================================\n");
    std::printf(" PHASIS -- Fase 1: nucleo em espaco plano\n");
    std::printf("=====================================================\n");

    // -----------------------------------------------------------------
    // T1: esfera homogenea, raio central.
    //     tau = N_A rho0 sigma * 2R
    // -----------------------------------------------------------------
    std::printf("\n  T1  esfera homogenea, b = 0\n");
    {
        const double rho0 = 3.7, R = 4.2e8, sig = 1.3e-33;
        const UniformBall ball(rho0, R);
        const PowerLawCrossSection xs(sig, 1.0, 0.0);   // alpha = 0: constante

        Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = 0.0;
        const Result r = trace_ray(ray, flat, ball, xs, opts_tight());

        check("tau", r.tau, N_A*rho0*sig*2.0*R, 1.0e-10);
        check("coluna X [g/cm^2]", r.column_density, rho0*2.0*R, 1.0e-10);
        check("comprimento proprio [cm]", r.path_length_cm, 2.0*R, 1.0e-10);
        check("r_min", r.r_min_cm, 0.0, 1.0e-14);
        check("P_surv", r.P_surv, std::exp(-N_A*rho0*sig*2.0*R), 1.0e-10);
    }

    // -----------------------------------------------------------------
    // T2: mesma esfera, b arbitrario.
    //     corda = 2 sqrt(R^2 - b^2)
    // -----------------------------------------------------------------
    std::printf("\n  T2  esfera homogenea, varredura em b/R\n");
    {
        const double rho0 = 3.7, R = 4.2e8, sig = 1.3e-33;
        const UniformBall ball(rho0, R);
        const PowerLawCrossSection xs(sig, 1.0, 0.0);

        for (double frac : {0.0, 0.1, 0.5, 0.9, 0.99}) {
            const double b = frac*R;
            Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = b;
            const Result r = trace_ray(ray, flat, ball, xs, opts_tight());

            const double corda = 2.0*std::sqrt(R*R - b*b);
            char nome[64];
            std::snprintf(nome, sizeof(nome), "tau  (b/R = %.2f)", frac);
            check(nome, r.tau, N_A*rho0*sig*corda, 1.0e-8);
            std::snprintf(nome, sizeof(nome), "r_min  (b/R = %.2f)", frac);
            check(nome, r.r_min_cm, b, 1.0e-14);
        }
    }

    // -----------------------------------------------------------------
    // T3: raio que passa por fora. tau = 0 exato, sem NaN.
    // -----------------------------------------------------------------
    std::printf("\n  T3  raio fora da materia (b > R)\n");
    {
        const double rho0 = 3.7, R = 4.2e8, sig = 1.3e-33;
        const UniformBall ball(rho0, R);
        const PowerLawCrossSection xs(sig, 1.0, 0.0);

        for (double frac : {1.0, 1.000001, 2.0, 1.0e6}) {
            Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = frac*R;
            const Result r = trace_ray(ray, flat, ball, xs, opts_tight());
            char nome[64];
            std::snprintf(nome, sizeof(nome), "tau  (b/R = %g)", frac);
            check(nome, r.tau, 0.0, 0.0);
            std::snprintf(nome, sizeof(nome), "P_surv  (b/R = %g)", frac);
            check(nome, r.P_surv, 1.0, 0.0);
            std::snprintf(nome, sizeof(nome), "finito  (b/R = %g)", frac);
            check_bool(nome, std::isfinite(r.tau), true);
            std::snprintf(nome, sizeof(nome), "crosses_matter (b/R = %g)", frac);
            check_bool(nome, r.crosses_matter, false);
        }
    }

    // -----------------------------------------------------------------
    // T4: halo em lei de potencia com p = 2.
    //
    // rho = rho0 (r/r0)^-2 entre r_in e r_out. Com dl = r dr/sqrt(r^2-b^2),
    //
    //   tau = 2 N_A sigma rho0 r0^2 * INT_{r_lo}^{r_out} dr/(r sqrt(r^2-b^2))
    //       = 2 N_A sigma rho0 r0^2 * (1/b) [ arccos(b/r_out) - arccos(b/r_lo) ]
    //
    // com r_lo = max(b, r_in).  No limite b -> 0 isto vira
    //
    //   tau = 2 N_A sigma rho0 r0^2 ( 1/r_in - 1/r_out ).
    // -----------------------------------------------------------------
    std::printf("\n  T4  halo p = 2 (forma fechada)\n");
    {
        const double rho0 = 2.5, r0 = 1.0e8, p = 2.0;
        const double r_in = 3.0e7, r_out = 9.0e8, sig = 7.0e-34;
        const PowerLawHalo halo(rho0, r0, p, r_in, r_out);
        const PowerLawCrossSection xs(sig, 1.0, 0.0);

        const double A = 2.0*N_A*sig*rho0*r0*r0;

        {   // b = 0
            Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = 0.0;
            const Result r = trace_ray(ray, flat, halo, xs, opts_tight());
            check("tau  (b = 0)", r.tau, A*(1.0/r_in - 1.0/r_out), 1.0e-9);
        }

        for (double b : {1.0e7, 3.0e7, 1.0e8, 5.0e8}) {
            Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = b;
            const Result r = trace_ray(ray, flat, halo, xs, opts_tight());

            const double r_lo = std::max(b, r_in);
            const double esperado =
                A*(1.0/b)*(std::acos(b/r_out) - std::acos(b/r_lo));

            char nome[64];
            std::snprintf(nome, sizeof(nome), "tau  (b = %.1e cm)", b);
            check(nome, r.tau, esperado, 1.0e-8);
        }
    }

    // -----------------------------------------------------------------
    // T5: sanity check de unidades com a Terra.
    // -----------------------------------------------------------------
    std::printf("\n  T5  Terra: coerencia de unidades\n");
    {
        const double R = units::earth_radius_cm;
        const double rho0 = units::earth_mean_density;
        const double sig = 1.0e-33;

        const UniformBall terra(rho0, R);
        const PowerLawCrossSection xs(sig, 1.0, 0.0);

        Ray ray; ray.E_inf_GeV = 1.0e6; ray.b_cm = 0.0;
        const Result r = trace_ray(ray, flat, terra, xs, opts_tight());

        std::printf("    R = %.4g cm,  rho = %.3g g/cm^3,  sigma = %.1e cm^2\n",
                    R, rho0, sig);
        check("coluna X [g/cm^2]  (esperado ~7e9)", r.column_density, 7.0208e9, 1.0e-4);
        check("tau                (esperado ~4)",   r.tau,            4.2281,   1.0e-4);
        check("P_surv",                             r.P_surv, std::exp(-r.tau), 1.0e-12);
    }

    // -----------------------------------------------------------------
    // T6: tabela de secao de choque, interpolacao log-log.
    // -----------------------------------------------------------------
    std::printf("\n  T6  TableCrossSection: log-log exato em lei de potencia\n");
    {
        const double s0 = 3.3e-35, E0 = 1.0e3, alpha = 0.402;
        const std::string caminho = "tests_tmp_sigma.csv";

        {   // tabela sintetica, grade grossa de proposito
            std::ofstream out(caminho);
            out << "# tabela sintetica: sigma = " << s0 << " * (E/" << E0
                << ")^" << alpha << "\n";
            out << "# E_GeV,sigma_cm2\n";
            out.precision(17);
            for (int i = 0; i <= 20; ++i) {
                const double E = E0*std::pow(10.0, i*11.0/20.0);
                out << E << "," << s0*std::pow(E/E0, alpha) << "\n";
            }
        }

        const TableCrossSection tab(caminho);
        const PowerLawCrossSection exato(s0, E0, alpha);

        check("numero de pontos", static_cast<double>(tab.size()), 21.0, 0.0);
        check_bool("cabecalho preservado", tab.header().size() == 2, true);

        // pontos INTERMEDIARIOS, longe dos nos da tabela
        double pior = 0.0;
        for (int i = 0; i < 97; ++i) {
            const double E = tab.E_min()*std::pow(tab.E_max()/tab.E_min(),
                                                  (i + 0.5)/97.0);
            const double a = tab.sigma_tot(E);
            const double b = exato.sigma_tot(E);
            pior = std::max(pior, std::fabs(a - b)/b);
        }
        check("pior erro em 97 pontos intermediarios", pior, 0.0, 1.0e-12);

        // fora da faixa deve LANCAR, nao extrapolar
        bool lancou_baixo = false, lancou_alto = false;
        try { tab.sigma_tot(0.5*tab.E_min()); } catch (const std::out_of_range&) { lancou_baixo = true; }
        try { tab.sigma_tot(2.0*tab.E_max()); } catch (const std::out_of_range&) { lancou_alto = true; }
        check_bool("lanca abaixo da faixa", lancou_baixo, true);
        check_bool("lanca acima da faixa",  lancou_alto,  true);

        std::remove(caminho.c_str());
    }

    // -----------------------------------------------------------------
    std::printf("\n=====================================================\n");
    std::printf("  %d verificacoes, %d falhas\n", total, falhas);
    std::printf("  RESULTADO: %s\n", falhas == 0 ? "OK" : "FALHA");
    std::printf("=====================================================\n\n");
    return falhas == 0 ? 0 : 1;
}
