// =====================================================================
// Corrente neutra de ponta a ponta -- T44 a T50.
//
// Fecha o P6: o dipolo produz sigma_NC(E) e dsigma_NC/dy(E,y), e o
// PHASIS consome os dois na cascata.
//
// O que mudou no PHASIS para isso ser possivel: CascadeKernel LANCAVA
// para qualquer secao de choque cuja forma em y dependesse de E, e
// deferia o caso para "Fase 4b". A tabela real DEPENDE -- medido aqui em
// T46 -- entao o caminho teve de ser construido.
// =====================================================================

#include "phasis/cascade.hpp"
#include "phasis/cross_section.hpp"
#include "phasis/density.hpp"
#include "phasis/metric.hpp"
#include "phasis/trace.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <memory>
#include <string>
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
    std::printf("    %-48s %-16.9g ref %-16.9g rel %8.2e  %s\n",
                n, g, r, rel, ok ? "OK" : "<-- FALHA");
}

void chk_bool(const char* n, bool g, bool r)
{
    ++total;
    const bool ok = (g == r);
    if (!ok) ++falhas;
    std::printf("    %-48s %-5s esperado %-5s              %s\n",
                n, g ? "sim" : "nao", r ? "sim" : "nao", ok ? "OK" : "<-- FALHA");
}

bool lanca(const std::function<void()>& f)
{
    try { f(); return false; } catch (const std::exception&) { return true; }
}

const char* kCC  = "dipole/data/sigma_nuN_CC_GBW.dat";
const char* kNC  = "dipole/data/sigma_nuN_NC_GBW.dat";
const char* kDY  = "dipole/data/dsigma_dy_NC_GBW.dat";

bool existe(const char* p) { std::ifstream f(p); return static_cast<bool>(f); }

void escreve_meta(const std::string& path, const std::string& projectile,
                  const std::string& target, const std::string& current)
{
    std::ofstream f(path);
    f << "# convention_y = (E_in - E_out)/E_in\n"
      << "# target       = " << target << "\n"
      << "# projectile   = " << projectile << "\n"
      << "# current      = " << current << "\n"
      << "# units_sigma  = cm^2\n"
      << "# units_E      = GeV\n"
      << "# M_Z_GeV      = 91.1876\n"
      << "# dipole_model = teste\n"
      << "# generated_by = tests/test_phase8.cpp\n";
    f.precision(17);
    for (int i = 0; i < 20; ++i) {
        const double E = 1.0e3*std::pow(1.0e10, i/19.0);
        f << E << " " << 1.0e-33*std::pow(E/1.0e6, 0.3) << "\n";
    }
}

} // namespace


int main()
{
    std::printf("\n=====================================================\n");
    std::printf("  CORRENTE NEUTRA DE PONTA A PONTA -- T44 a T50\n");
    std::printf("=====================================================\n");

    if (!existe(kCC) || !existe(kNC) || !existe(kDY)) {
        std::printf("\n  <-- FALHA: faltam tabelas do dipolo.\n");
        std::printf("      Gere com, dentro de dipole/:\n");
        std::printf("        ./build/sigma_nuN --model GBW --current CC --NE 300 ...\n");
        std::printf("        ./build/sigma_nuN --model GBW --current NC --NE 300 --nY 120 ...\n");
        return 1;
    }

    // =================================================================
    std::printf("\n  T45  COMPOSICAO: somar CC com NC exige o guarda OPOSTO ao\n");
    std::printf("       de comparar. assert_comparable pede current IGUAL.\n");
    {
        const std::string a = "/tmp/phasis_t45_cc.dat";
        const std::string b = "/tmp/phasis_t45_nc.dat";

        escreve_meta(a, "nu", "isoscalar_nucleon", "CC");
        escreve_meta(b, "nu", "isoscalar_nucleon", "NC");
        chk_bool("CC + NC: passa",
                 lanca([&]{ assert_composable(TableCrossSection(a), TableCrossSection(b)); }), false);

        escreve_meta(b, "nu", "isoscalar_nucleon", "CC");
        chk_bool("CC + CC: lanca (contaria duas vezes)",
                 lanca([&]{ assert_composable(TableCrossSection(a), TableCrossSection(b)); }), true);

        escreve_meta(b, "nubar", "isoscalar_nucleon", "NC");
        chk_bool("nu + nubar: lanca",
                 lanca([&]{ assert_composable(TableCrossSection(a), TableCrossSection(b)); }), true);

        // E o guarda de comparar continua rejeitando CC contra NC.
        escreve_meta(b, "nu", "isoscalar_nucleon", "NC");
        chk_bool("assert_comparable(CC, NC): lanca",
                 lanca([&]{ assert_comparable(TableCrossSection(a), TableCrossSection(b)); }), true);

        // As tabelas de producao.
        chk_bool("producao CC + NC: passa",
                 lanca([&]{ assert_composable(TableCrossSection(kCC), TableCrossSection(kNC)); }), false);
    }

    const TableCrossSection cc_tab(kCC);
    const TableCrossSection nc_tab(kNC);

    // =================================================================
    std::printf("\n  T44  sigma_NC / sigma_CC nas tabelas reais\n");
    {
        // O valor previsto sai das somas efetivas de carga, e nao depende
        // de nenhum acoplamento:
        //
        //   NC: Soma_q (g_V^2+g_A^2) sobre u,d,s,c = 1.3127   (o C de KK)
        //   CC: 2 canais x (g_V^2+g_A^2) = 2         -> 4
        //
        // razao = (1.3127/4) (M_Z/M_W)^2 = 0.4224, que e o valor do SM.
        const double s2w = 0.23122;
        const double gVu = 0.5 - (4.0/3.0)*s2w, gAu = 0.5;
        const double gVd = -0.5 + (2.0/3.0)*s2w, gAd = -0.5;
        const double S = 2.0*(gVu*gVu + gAu*gAu) + 2.0*(gVd*gVd + gAd*gAd);
        const double previsto = (S/4.0)*std::pow(91.1876/80.379, 2);
        std::printf("    soma efetiva NC = %.5f (u,d,s,c);  CC = 4\n", S);
        std::printf("    previsto (S/4)(M_Z/M_W)^2 = %.5f ; SM conhecido ~0.42\n", previsto);

        std::printf("    %-13s %-15s %-15s %s\n", "E_GeV", "sigma_CC", "sigma_NC", "razao");
        double r_alto = 0.0; int n_alto = 0;
        bool monotona = true, anterior_ok = false; double r_ant = 0.0;
        for (double E : {1.0e3, 1.0e5, 1.0e7, 1.0e9, 1.0e11, 1.0e13}) {
            const double a = cc_tab.sigma_tot(E), b = nc_tab.sigma_tot(E);
            std::printf("    %-13.3g %-15.5e %-15.5e %.5f\n", E, a, b, b/a);
            if (E >= 1.0e9) { r_alto += b/a; ++n_alto; }
            if (anterior_ok && !(b/a >= r_ant - 1.0e-12)) monotona = false;
            r_ant = b/a; anterior_ok = true;
        }
        r_alto /= n_alto;
        // A razao SOBE com a energia porque o xF3 do CC domina em baixo
        // (e de valencia) e some em cima. O valor assintotico e o que a
        // contagem de cargas preve.
        chk("razao assintotica (media em E >= 1e9)", r_alto, previsto, 0.05);
        chk_bool("razao cresce monotonicamente com E", monotona, true);
    }

    // =================================================================
    std::printf("\n  T46  A FORMA EM y DEPENDE DE E -- e por isso o kernel\n");
    std::printf("       de forma fixa nao servia\n");

    auto cc_ptr = std::make_shared<TableCrossSection>(kCC);
    TableDifferentialCrossSection::Expect esp;   // convention_y e current = NC
    const TableDifferentialCrossSection nc_dif(kDY, esp, cc_ptr);
    {
        chk_bool("o leitor recusa forma E-independente", nc_dif.shape_is_E_independent(), false);
        std::printf("    g(y) = (dsigma/dy)/sigma_nc, adimensional:\n");
        std::printf("    %-10s", "y");
        for (double E : {1.0e5, 1.0e8, 1.0e11}) std::printf("  E=%-11.0e", E);
        std::printf("\n");
        double maior_drift = 0.0;
        for (double y : {1.0e-3, 1.0e-2, 1.0e-1, 0.5}) {
            std::printf("    %-10.0e", y);
            double lo = 1.0e300, hi = 0.0;
            for (double E : {1.0e5, 1.0e8, 1.0e11}) {
                const double g = nc_dif.dsigma_nc_dy(E, y)/nc_dif.sigma_nc(E);
                std::printf("  %-13.5f", g);
                lo = std::min(lo, g); hi = std::max(hi, g);
            }
            std::printf("\n");
            if (lo > 0.0) maior_drift = std::max(maior_drift, hi/lo);
        }
        std::printf("    maior variacao da forma entre 1e5 e 1e11: fator %.1f\n", maior_drift);
        ++total;
        if (maior_drift > 2.0) {
            std::printf("    OK: a forma varia por fator %.0f -- congela-la seria erro\n", maior_drift);
        } else { ++falhas;
            std::printf("    <-- FALHA: a forma quase nao varia; revisar a premissa\n"); }

        // Resolucao da grade em y: INT dy contra a tabela total.
        std::printf("    resolucao da grade em y (INT dy vs tabela total de sigma_NC):\n");
        double pior = 0.0;
        for (double E : {1.0e4, 1.0e6, 1.0e9, 1.0e12}) {
            const double r = nc_dif.sigma_nc(E)/nc_tab.sigma_tot(E);
            std::printf("      E=%-11.3g  %.5f\n", E, r);
            pior = std::max(pior, std::fabs(r - 1.0));
        }
        ++total;
        if (pior <= 0.02) {
            std::printf("    OK: o trapezio sobre a grade recupera sigma em %.2f%%\n", 100*pior);
        } else { ++falhas;
            std::printf("    <-- FALHA: desvio de %.2f%% -- grade em y grossa demais\n", 100*pior); }
    }

    // =================================================================
    std::printf("\n  T47  IDENTIDADE DE CONSISTENCIA com a tabela real\n");
    const EnergyGrid grid(1.0e7, 1.0e11, 24);
    const CascadeKernel ker(grid, nc_dif, 1.0e-12, 4096, /*n_E_nodes=*/6);
    {
        chk_bool("o kernel sabe que a forma depende de E", ker.shape_depends_on_E(), true);
        std::printf("    %d nos em E, faixa [%.3g, %.3g] GeV\n",
                    ker.n_E_nodes(), ker.kernel_E_min(), ker.kernel_E_max());
        chk("pior residuo nos nos de E", ker.worst_residual(), 0.0, 1.0e-14);

        // O ponto da interpolacao LINEAR: a identidade e linear, logo
        // sobrevive entre os nos, nao so neles. Testado em E que nao
        // coincide com no nenhum.
        double pior_entre = 0.0;
        for (double t : {0.13, 0.37, 0.61, 0.88}) {
            const double E = ker.kernel_E_min()
                *std::pow(ker.kernel_E_max()/ker.kernel_E_min(), t);
            for (int j = 0; j < grid.n(); ++j) {
                double soma = ker.G_leak(j, E);
                for (int i = 0; i <= j; ++i) soma += ker.G(i, j, E);
                pior_entre = std::max(pior_entre, std::fabs(soma - 1.0));
            }
        }
        chk("identidade ENTRE os nos (interpolacao linear)", pior_entre, 0.0, 1.0e-14);
        std::printf("    desvio de norma (trapezio da tabela vs adaptativo): %.3f%%\n",
                    100*ker.worst_norm_mismatch());

        // G(i,j) sem energia tem de LANCAR quando a forma depende de E:
        // devolver o primeiro no seria um numero plausivel e errado.
        chk_bool("G(i,j) sem E lanca", lanca([&]{ (void)ker.G(0, 5); }), true);
        chk_bool("G_leak(j) sem E lanca", lanca([&]{ (void)ker.G_leak(5); }), true);
        chk_bool("fora da faixa em E lanca",
                 lanca([&]{ (void)ker.G(0, 5, 1.0e-3*ker.kernel_E_min()); }), true);
    }

    // =================================================================
    std::printf("\n  T48  CONVERGENCIA nos nos de E\n");
    {
        const CascadeKernel k12(grid, nc_dif, 1.0e-12, 4096, 12);
        double pior = 0.0;
        for (double t : {0.2, 0.5, 0.8}) {
            const double E = ker.kernel_E_min()
                *std::pow(ker.kernel_E_max()/ker.kernel_E_min(), t);
            for (int j = 0; j < grid.n(); ++j)
                for (int i = 0; i <= j; ++i) {
                    const double a = ker.G(i, j, E), b = k12.G(i, j, E);
                    if (b > 1.0e-6) pior = std::max(pior, std::fabs(a/b - 1.0));
                }
        }
        std::printf("    6 nos contra 12 nos, maior diferenca relativa em G: %.3e\n", pior);
        ++total;
        if (pior <= 0.05) std::printf("    OK: 6 nos ja resolvem a forma em ln E\n");
        else { ++falhas; std::printf("    <-- FALHA: precisa de mais nos em E\n"); }
    }

    // =================================================================
    std::printf("\n  T49  os kernels ANALITICOS continuam pelo caminho antigo\n");
    {
        ConstantY xs(1.0e-33, 5.0e-34, 1.0e6, 0.36, 0.36);
        const CascadeKernel k(grid, xs, 1.0e-12);
        chk_bool("forma E-independente: um no so", k.shape_depends_on_E(), false);
        chk_bool("G(i,j) sem E funciona", lanca([&]{ (void)k.G(0, 3); }), false);
        chk("residuo", k.worst_residual(), 0.0, 1.0e-14);
        // Z analitico: g(y) = 1 => Z = 1/gamma.
        //
        // Aqui NAO cabe tolerancia fixa: Z_discreto e uma soma de
        // Riemann sobre a grade, e a 6 bins/decada ele erra ~2% por
        // discretizacao, nao por bug. O que vale testar e a ORDEM: o
        // erro tem de cair como o quadrado do espacamento.
        std::printf("    convergencia de Z_discrete(2) para 1/gamma = 0.5:\n");
        double ant = 0.0, pior_ordem = 1.0e9;
        for (int bpd : {6, 12, 24, 48}) {
            const EnergyGrid gz(1.0e7, 1.0e11, 4*bpd);
            const CascadeKernel kz(gz, xs, 1.0e-12);
            const double err = std::fabs(kz.Z_discrete(2.0)/0.5 - 1.0);
            std::printf("      %3d bins/decada  Z = %.8f  erro %.3e", bpd,
                        kz.Z_discrete(2.0), err);
            if (ant > 0.0) {
                const double ordem = std::log2(ant/err);
                pior_ordem = std::min(pior_ordem, ordem);
                std::printf("  ordem %.2f", ordem);
            }
            std::printf("\n");
            ant = err;
        }
        ++total;
        if (pior_ordem >= 1.8) {
            std::printf("    OK: convergencia de segunda ordem (pior %.2f)\n", pior_ordem);
        } else { ++falhas;
            std::printf("    <-- FALHA: ordem %.2f, abaixo de 2\n", pior_ordem); }
    }

    // =================================================================
    std::printf("\n  T50  TRANSPORTE de ponta a ponta com as tabelas reais\n");
    {
        const double r_s = 1.0e6;
        const Schwarzschild schw(r_s);
        const PowerLawHalo halo(1.0e-2, 1.0e8, 2.0, 3.0*r_s, 1.0e4*r_s);
        IntegratorOpts o; o.ode_rel_tol = 1.0e-12;

        Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = 5.0*r_s;

        std::vector<double> phi0(grid.n());
        for (int i = 0; i < grid.n(); ++i) phi0[i] = grid.powerlaw_bin_integral(i, 2.0);
        double s0 = 0.0; for (double v : phi0) s0 += v;

        // (a) SEM absorcao CC: a NC so redistribui e vaza, entao
        //     soma(phi) + vazamento tem de conservar EXATAMENTE.
        auto cc_zero = std::make_shared<ScaledCrossSection>(cc_ptr, 0.0);
        const TableDifferentialCrossSection nc_sem_cc(kDY, esp, cc_zero);
        const CascadeKernel ker0(grid, nc_sem_cc, 1.0e-12, 4096, 6);
        const CascadeResult a = transport_cascade(phi0, ker0, ray, schw, halo, o);
        double sa = 0.0; for (double v : a.phi) sa += v;
        std::printf("    sem CC: soma(phi) = %.10f, vazamento = %.10f, phi0 = %.10f\n",
                    sa, a.leakage_total, s0);
        chk("conservacao de numero sem CC", (sa + a.leakage_total)/s0, 1.0, 1.0e-9);

        // (b) COM absorcao CC, contra a mesma geometria sem regeneracao.
        const CascadeResult b = transport_cascade(phi0, ker, ray, schw, halo, o);
        double sb = 0.0; for (double v : b.phi) sb += v;

        ZeroNC zero(1.0, 1.0e6, 0.0);   // sem regeneracao, sigma_cc irrelevante
        (void)zero;
        // referencia: mesma absorcao, sem ganho -- exp(-tau) bin a bin
        const SumCrossSection* dummy = nullptr; (void)dummy;
        SumCrossSection tot;
        tot.add(cc_ptr);
        tot.add(std::make_shared<TableCrossSection>(kNC));
        double s_ref = 0.0;
        for (int i = 0; i < grid.n(); ++i) {
            Ray r2 = ray; r2.E_inf_GeV = grid.center(i);
            const Result t = trace_ray(r2, schw, halo, tot, o);
            s_ref += phi0[i]*std::exp(-t.tau);
        }
        std::printf("    com CC: soma(phi) = %.6e ; so absorcao = %.6e ; ganho = %+.2f%%\n",
                    sb, s_ref, 100*(sb/s_ref - 1.0));
        chk_bool("regeneracao NC AUMENTA o fluxo sobrevivente", sb > s_ref, true);
        chk_bool("o vazamento e rastreado e positivo", b.leakage_total > 0.0, true);
        chk_bool("tolerancia atendida", b.tolerance_met, true);
    }

    std::printf("\n=====================================================\n");
    std::printf("  %d verificacoes, %d falhas\n", total, falhas);
    std::printf("  RESULTADO: %s\n", falhas == 0 ? "OK" : "FALHA");
    std::printf("=====================================================\n\n");
    return falhas == 0 ? 0 : 1;
}
