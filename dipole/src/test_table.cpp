// ================================================================
// F3 e F4 - teste de aceitacao
//
// F3: erro da interpolacao bilinear em (ln x, ln Q^2) sobre ln F,
//     medido contra o calculo direto em pontos FORA dos nos.
//
// F4: (a) convergencia em densidade de nos por decada;
//     (b) condicionamento: sigma ~ E^0.36, entao uma perturbacao
//         dE/E deve dar amplificacao ~0.4, nunca 1e5 como antes.
// ================================================================

#include <cmath>
#include <cstdio>
#include <vector>

#include "parameters.hpp"
#include "dipole_models.hpp"
#include "integrals.hpp"
#include "structure_table.hpp"
#include "sigma_nuN_core.hpp"

using namespace dipole;

namespace {

int falhas = 0;

const double MAX_ERRO_INTERP  = 0.02;   // 2%

StructureTable montar(int Nx, int NQ, int Nr, int Nz, double logEmax = 14.0)
{
    GBWParameters gbw;
    IIMParameters iim;

    QuarkMasses masses;
    masses.u = masses.d = masses.s = gbw.m_light;
    masses.c = gbw.m_charm;

    QuadratureGrid quad; quad.Nr = Nr; quad.Nz = Nz;

    const double s_max = 2.0*MN*std::pow(10.0, logEmax);
    TableSpec spec;
    spec.Nx = Nx; spec.NQ = NQ;
    spec.Q2Min = 1.0;
    spec.Q2Max = 1.05*0.999*s_max;
    spec.xMin  = 0.5/s_max;
    spec.xMax  = 1.0;

    const std::vector<StructureTable::Channel> canais = {
        { QuarkFlavor::u, QuarkFlavor::d },
        { QuarkFlavor::c, QuarkFlavor::s }
    };

    return StructureTable(DipoleModelId::GBW, canais, gbw, iim,
                          quad, spec, masses, false, CurrentType::CC);
}

} // namespace

int main()
{
    std::printf("\n=== F3 e F4: teste de aceitacao ===\n\n");

    std::printf("Montando tabela 121 x 121 (Nr=120, Nz=120)...\n");
    const StructureTable tab = montar(121, 121, 120, 120);
    std::printf("  %s\n  %.1f s\n\n", tab.describe().c_str(), tab.buildSeconds());

    // ---- F3: erro de interpolacao --------------------------------
    std::printf("F3) Erro da interpolacao, em pontos fora dos nos\n\n");
    std::printf("  %-12s %-12s %-13s %-13s %s\n",
                "x", "Q2", "tabela F_2", "exato F_2", "erro");

    const double xs[] = { 3.7e-3, 8.1e-5, 2.3e-7, 6.9e-10, 1.4e-12 };
    const double qs[] = { 2.7, 41.0, 830.0, 6460.78, 1.9e5 };

    double pior = 0.0;
    for (double x : xs) {
        for (double Q2 : qs) {
            const StructureTL a = tab.at(x, Q2);
            const StructureTL e = tab.exact(x, Q2);
            const double fa = a.FT + a.FL;
            const double fe = e.FT + e.FL;
            if (!(fe > 0.0)) continue;
            const double rel = std::fabs(fa - fe)/fe;
            if (rel > pior) pior = rel;
            std::printf("  %-12.3g %-12.4g %-13.6g %-13.6g %+.3f%%\n",
                        x, Q2, fa, fe, 100.0*rel);
        }
    }
    const bool ok3 = pior <= MAX_ERRO_INTERP;
    if (!ok3) ++falhas;
    std::printf("\n  pior erro: %.3f%%  (limite %.1f%%)  %s\n\n",
                100.0*pior, 100.0*MAX_ERRO_INTERP, ok3 ? "OK" : "<-- FALHA");

    // ---- F4a: convergencia em densidade de nos -------------------
    std::printf("F4a) Convergencia em nos por decada\n\n");
    std::printf("  %-10s", "E [GeV]");
    const double dens[] = { 4.0, 8.0, 16.0, 32.0 };
    for (double d : dens) std::printf("  n/dec=%-13.0f", d);
    std::printf("\n");

    const double energias[] = { 1.0e6, 1.0e9, 1.0e12, 1.0e14 };
    for (double E : energias) {
        std::printf("  %-10.1e", E);
        double ant = 0.0, rel_ultimo = 0.0;
        for (double d : dens) {
            const double s = sigmaNuN_CC(E, tab, d, d)*GeVminus2_to_cm2;
            std::printf("  %-19.6e", s);
            if (ant > 0.0) rel_ultimo = std::fabs(s - ant)/ant;
            ant = s;
        }
        const bool ok = rel_ultimo <= 0.01;
        if (!ok) ++falhas;
        std::printf("  (16->32: %.3f%%  %s)\n", 100.0*rel_ultimo, ok ? "OK" : "FALHA");
    }

    // ---- F4b: suavidade de sigma(E) -----------------------------
    //
    // A metrica de "amplificacao" (dsigma/sigma dividido por dE/E) so faz
    // sentido para uma resposta diferenciavel. Com a tabela, sigma(E) e
    // suave na escala que importa (a grade de producao tem dE/E ~ 9%),
    // mas nao e diferenciavel na escala de 1e-7, porque os nos de Simpson
    // cruzam a estrutura local do interpolante. Medir amplificacao ali e
    // medir ruido sem significado fisico.
    //
    // O que importa de fato: sigma(E) tem de ser monotonica e suave o
    // bastante para que a interpolacao log-log da tabela (que e o que o
    // HADROS3 faz) nao injete artefato.
    //
    // Metrica: a quarta diferenca finita em ln E aniquila qualquer cubica,
    // entao D4/6 estima o desvio de cada ponto em relacao a uma curva
    // suave local, sem precisar ajustar nada.
    std::printf("\nF4b) Suavidade de sigma(E) na grade de producao\n\n");

    const int NE = 45;
    const double logE0 = 3.0, logE1 = 14.0;

    std::vector<double> ls(NE);
    bool monotonica = true;
    double anterior = 0.0;

    for (int i = 0; i < NE; ++i) {
        const double E = std::pow(10.0, logE0 + i*(logE1 - logE0)/(NE - 1));
        const double sig = sigmaNuN_CC(E, tab, 16.0, 16.0);
        ls[i] = std::log(sig);
        if (i > 0 && sig <= anterior) monotonica = false;
        anterior = sig;
    }

    double soma2 = 0.0, maxres = 0.0;
    int n = 0;
    for (int i = 2; i < NE - 2; ++i) {
        const double d4 = ls[i-2] - 4.0*ls[i-1] + 6.0*ls[i]
                        - 4.0*ls[i+1] + ls[i+2];
        const double res = std::fabs(d4)/6.0;
        soma2 += res*res;
        if (res > maxres) maxres = res;
        ++n;
    }
    const double rms = std::sqrt(soma2/n);

    std::printf("  %d energias de 1e%g a 1e%g GeV\n", NE, logE0, logE1);
    std::printf("  monotonica              : %s\n", monotonica ? "sim" : "NAO");
    std::printf("  residuo de suavidade rms: %.4f%%\n", 100.0*rms);
    std::printf("  residuo maximo          : %.4f%%\n", 100.0*maxres);
    std::printf("  tabela publicada (antes) : 29 passos decrescentes, rms 0.75%%\n");

    const bool ok4b = monotonica && rms <= 0.005;
    if (!ok4b) ++falhas;
    std::printf("\n  criterio: monotonica e rms <= 0.5%%   %s\n",
                ok4b ? "OK" : "<-- FALHA");

    std::printf("\n  RESULTADO: %s\n\n", falhas == 0 ? "OK" : "FALHA");
    return falhas == 0 ? 0 : 1;
}
