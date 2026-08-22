// =====================================================================
// Fase 5 -- infraestrutura de varredura. T27, T28, T29.
// (O ensemble de emissao ainda nao esta decidido; nada aqui depende dele.)
// =====================================================================

#include "phasis/cross_section.hpp"
#include "phasis/density.hpp"
#include "phasis/metric.hpp"
#include "phasis/shell.hpp"
#include "phasis/trace.hpp"
#include "phasis/units.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace phasis;

namespace {

int falhas = 0, total = 0;

void chk(const char* n, double g, double r, double t)
{
    ++total;
    const double d = (std::fabs(r) > 0.0) ? std::fabs(r) : 1.0;
    const double rel = std::fabs(g - r)/d;
    const bool ok = (rel <= t) && std::isfinite(g);
    if (!ok) ++falhas;
    std::printf("    %-46s %-18.11g ref %-18.11g rel %8.2e  %s\n",
                n, g, r, rel, ok ? "OK" : "<-- FALHA");
}

void chk_bool(const char* n, bool g, bool r)
{
    ++total;
    const bool ok = (g == r);
    if (!ok) ++falhas;
    std::printf("    %-46s %-18s ref %-18s                %s\n",
                n, g ? "true" : "false", r ? "true" : "false", ok ? "OK" : "<-- FALHA");
}

} // namespace

int main()
{
    const Minkowski flat;
    const double NA = units::N_A;

    std::printf("\n=====================================================\n");
    std::printf(" PHASIS -- Fase 5: localizacao da casca tau = 1\n");
    std::printf("=====================================================\n");

    // =================================================================
    std::printf("\n  T27  casca analitica: Minkowski + esfera + sigma constante\n");
    {
        // tau(b) = 2 N_A sigma rho sqrt(R^2 - b^2)
        //   =>  b* = sqrt( R^2 - (1/(2 N_A sigma rho))^2 )
        const double R = 4.0e8, rho0 = 2.0, sig = 3.0e-33;
        const UniformBall bola(rho0, R);
        const PowerLawCrossSection xs(sig, 1.0, 0.0);

        auto tau_de_b = [&](double b) {
            Ray r; r.E_inf_GeV = 1.0e9; r.b_cm = b;
            IntegratorOpts o; o.rel_tol = 1.0e-13;
            return trace_ray(r, flat, bola, xs, o).tau;
        };

        const double c = 1.0/(2.0*NA*sig*rho0);
        const double b_exato = std::sqrt(R*R - c*c);

        ShellOpts so; so.n_coarse = 64;
        const ShellResult sr = locate_shell(tau_de_b, 0.0, R, so);
        std::printf("    %s\n", sr.diagnostico().c_str());

        chk_bool("casca encontrada", sr.shell_found, true);
        chk_bool("raiz unica", !sr.multi_root, true);
        chk("b* localizado", sr.roots.empty() ? 0.0 : sr.roots[0], b_exato, 1.0e-12);
        chk("tau em b*", tau_de_b(sr.roots.empty() ? 0.0 : sr.roots[0]), 1.0, 1.0e-10);
    }

    // =================================================================
    std::printf("\n  T28  nao-monotonicidade: perfil OCO da duas raizes\n");
    {
        // Este teste existe porque bissecao cega e o unico modo de falha
        // que nao produz erro visivel: ela devolveria UMA raiz e nada
        // avisaria.
        //
        // Basta uma casca oca [r1, r2]. A corda vale
        //   b < r1 : 2[sqrt(r2^2-b^2) - sqrt(r1^2-b^2)]
        //   b > r1 : 2 sqrt(r2^2-b^2)
        // e CRESCE de 2(r2-r1) em b=0 ate 2 sqrt(r2^2-r1^2) em b=r1,
        // so entao caindo a zero. Com o alvo entre esses dois valores,
        // ha exatamente duas raizes.
        const double r1 = 1.0e8, r2 = 2.0e8;
        const double sig = 1.0e-33;
        const double corda0 = 2.0*(r2 - r1);                    // 2.0e8
        const double corda1 = 2.0*std::sqrt(r2*r2 - r1*r1);     // 3.464e8
        const double corda_alvo = 2.5e8;                        // entre as duas
        const double rho0 = 1.0/(NA*sig*corda_alvo);

        const PowerLawHalo oca(rho0, 1.0e8, 0.0, r1, r2);       // p = 0: densidade constante
        const PowerLawCrossSection xs(sig, 1.0, 0.0);

        auto tau_de_b = [&](double b) {
            Ray r; r.E_inf_GeV = 1.0e9; r.b_cm = b;
            IntegratorOpts o; o.rel_tol = 1.0e-13;
            return trace_ray(r, flat, oca, xs, o).tau;
        };

        std::printf("    corda: b=0 -> %.4e ;  b=r1 -> %.4e ;  alvo -> %.4e\n",
                    corda0, corda1, corda_alvo);
        std::printf("    tau:   b=0 -> %.5f ;  b=r1 -> %.5f\n",
                    tau_de_b(0.0), tau_de_b(r1*0.999999));

        ShellOpts so; so.n_coarse = 64;
        const ShellResult sr = locate_shell(tau_de_b, 0.0, r2, so);
        std::printf("    %s\n", sr.diagnostico().c_str());

        chk_bool("casca encontrada", sr.shell_found, true);
        chk_bool("MULTI_ROOT sinalizado", sr.multi_root, true);
        chk("numero de raizes", static_cast<double>(sr.roots.size()), 2.0, 0.0);

        // as duas contra a forma fechada
        if (sr.roots.size() == 2) {
            // raiz 1, com b < r1: 2[sqrt(r2^2-b^2) - sqrt(r1^2-b^2)] = corda_alvo
            double lo = 0.0, hi = r1;
            for (int i = 0; i < 200; ++i) {
                const double m = 0.5*(lo + hi);
                const double c = 2.0*(std::sqrt(r2*r2-m*m) - std::sqrt(r1*r1-m*m));
                // A corda CRESCE com b nesta faixa (de 2(r2-r1) ate
                // 2 sqrt(r2^2-r1^2)), entao c < alvo significa ir para b
                // MAIOR. Errei o sentido na primeira versao; o
                // localizador estava certo e a referencia e que estava
                // invertida -- tau na raiz que ele achou da 1 a 1e-13.
                if (c < corda_alvo) lo = m; else hi = m;
            }
            const double b1 = 0.5*(lo + hi);
            // raiz 2, com b > r1: 2 sqrt(r2^2-b^2) = corda_alvo
            const double b2 = std::sqrt(r2*r2 - 0.25*corda_alvo*corda_alvo);

            chk("raiz 1 vs forma fechada", sr.roots[0], b1, 1.0e-9);
            chk("raiz 2 vs forma fechada", sr.roots[1], b2, 1.0e-12);
            chk("tau na raiz 1", tau_de_b(sr.roots[0]), 1.0, 1.0e-9);
            chk("tau na raiz 2", tau_de_b(sr.roots[1]), 1.0, 1.0e-10);
        }

        // Contraste: bissecao cega no intervalo inteiro acharia UMA e
        // acharia a errada, sem sinal nenhum de que errou.
        std::printf("    (bissecao cega em [0, r2] devolveria so uma das duas,\n");
        std::printf("     sem nada indicando que a outra existe)\n");
    }

    // =================================================================
    std::printf("\n  T29  sem bracket: transparente e opaco em toda parte\n");
    {
        const double R = 4.0e8;
        const PowerLawCrossSection xs(1.0e-33, 1.0, 0.0);

        for (int caso = 0; caso < 2; ++caso) {
            const double rho0 = caso ? 1.0e3 : 1.0e-9;   // opaco / transparente
            const UniformBall bola(rho0, R);

            auto tau_de_b = [&](double b) {
                Ray r; r.E_inf_GeV = 1.0e9; r.b_cm = b;
                IntegratorOpts o; o.rel_tol = 1.0e-11;
                return trace_ray(r, flat, bola, xs, o).tau;
            };

            // faixa que NAO chega a b = R, onde tau cairia a zero
            ShellOpts so; so.n_coarse = 48;
            const ShellResult sr = locate_shell(tau_de_b, 0.0, 0.9*R, so);

            std::printf("    %-14s %s\n", caso ? "denso:" : "rarefeito:",
                        sr.diagnostico().c_str());
            chk_bool(caso ? "denso: shell_found = false"
                          : "rarefeito: shell_found = false", sr.shell_found, false);
            chk_bool(caso ? "denso: sem NaN" : "rarefeito: sem NaN",
                     std::isfinite(sr.tau_min) && std::isfinite(sr.tau_max), true);
            chk_bool(caso ? "denso: nenhuma raiz" : "rarefeito: nenhuma raiz",
                     sr.roots.empty(), true);
        }
    }

    std::printf("\n=====================================================\n");
    std::printf("  %d verificacoes, %d falhas\n", total, falhas);
    std::printf("  RESULTADO: %s\n", falhas == 0 ? "OK" : "FALHA");
    std::printf("=====================================================\n\n");
    return falhas == 0 ? 0 : 1;
}
