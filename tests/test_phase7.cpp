// =====================================================================
// Acoplamento GR x saturacao -- T36 a T39.
//
// A grandeza medida e propriedade so de sigma_tot. Nao depende de NC,
// nem de dsigma/dy, nem do acoplamento NC do dipole.
//
// A ORDEM IMPORTA: T39 roda cedo. Se R nao for invariante sob reescala
// de sigma, a medida esta contaminada pela normalizacao e comparar
// dipolo com colinear nao significa nada -- e ai nenhum dos outros
// numeros vale a pena olhar.
// =====================================================================

#include "phasis/cross_section.hpp"
#include "phasis/density.hpp"
#include "phasis/metric.hpp"
#include "phasis/trace.hpp"

#include <algorithm>
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
    std::printf("    %-46s %-18.12g ref %-18.12g rel %8.2e  %s\n",
                n, g, r, rel, ok ? "OK" : "<-- FALHA");
}

void chk_bool(const char* n, bool g, bool r)
{
    ++total;
    const bool ok = (g == r);
    if (!ok) ++falhas;
    std::printf("    %-46s %-5s esperado %-5s              %s\n",
                n, g ? "sim" : "nao", r ? "sim" : "nao", ok ? "OK" : "<-- FALHA");
}

// ---------------------------------------------------------------------
// R(b, E) = (tau_full - tau_frozen)/tau_frozen, com a MESMA geometria.
struct Rmed { double R, tau_frozen, tau_full, r_turn; bool ok; };

Rmed medir_R(const CrossSection& xs, const Metric& met,
             const DensityProfile& prof, double b, double E,
             const IntegratorOpts& o)
{
    Ray ray; ray.E_inf_GeV = E; ray.b_cm = b;
    const Result a = trace_ray(ray, met, prof, xs, o, false);
    const Result c = trace_ray(ray, met, prof, xs, o, true);
    Rmed m;
    m.tau_full   = a.tau;
    m.tau_frozen = c.tau;
    m.R = (c.tau > 0.0) ? (a.tau - c.tau)/c.tau : 0.0;
    m.r_turn = a.r_min_cm;
    m.ok = a.tolerance_met && c.tolerance_met && !a.captured;
    return m;
}

const char* kTabelaDipolo = "dipole/data/sigma_nuN_CC_GBW.dat";

// Tabela sintetica ONDULADA. Nao e uma lei de potencia: se fosse, o
// interpolador log-log a reproduziria exatamente e T39 nao exercitaria
// nada da interpolacao. As ondulacoes garantem que sigma(E) usado pela
// quadratura seja mesmo o resultado de interpolar.
std::string escreve_tabela_ondulada(const std::string& path)
{
    std::ofstream f(path);
    f << "# convention_y = (E_in - E_out)/E_in\n"
      << "# target       = isoscalar_nucleon\n"
      << "# projectile   = nu\n"
      << "# current      = CC\n"
      << "# units_sigma  = cm^2\n"
      << "# units_E      = GeV\n"
      << "# M_Z_GeV      = 91.1876\n"
      << "# dipole_model = sintetica-ondulada (teste)\n"
      << "# generated_by = tests/test_phase7.cpp\n";
    f.precision(17);
    const int n = 220;
    for (int i = 0; i < n; ++i) {
        const double E = 1.0e3*std::pow(1.0e11, double(i)/(n-1));
        const double L = std::log(E/1.0e6);
        const double s = 1.0e-33*std::pow(E/1.0e6, 0.30)
                         *(1.0 + 0.05*std::sin(1.7*L) + 0.02*std::cos(4.3*L));
        f << E << " " << s << "\n";
    }
    return path;
}

// Escreve uma tabela minima com metadados controlados, para exercitar
// a validacao cruzada.
void escreve_com_meta(const std::string& path, const std::string& projectile,
                      const std::string& target, const std::string& current,
                      bool completa = true)
{
    std::ofstream f(path);
    f << "# convention_y = (E_in - E_out)/E_in\n"
      << "# target       = " << target << "\n"
      << "# projectile   = " << projectile << "\n"
      << "# current      = " << current << "\n"
      << "# units_sigma  = cm^2\n"
      << "# units_E      = GeV\n";
    if (completa) {
        f << "# M_Z_GeV      = 91.1876\n"
          << "# dipole_model = teste\n"
          << "# generated_by = tests/test_phase7.cpp\n";
    }
    f.precision(17);
    for (int i = 0; i < 20; ++i) {
        const double E = 1.0e3*std::pow(1.0e10, i/19.0);
        f << E << " " << 1.0e-33*std::pow(E/1.0e6, 0.3) << "\n";
    }
}

bool lanca(const std::function<void()>& f)
{
    try { f(); return false; } catch (const std::exception&) { return true; }
}

} // namespace


int main()
{
    std::printf("\n=====================================================\n");
    std::printf("  ACOPLAMENTO GR x SATURACAO -- T36 a T39\n");
    std::printf("=====================================================\n");

    const double r_s = 1.0e6;
    const Schwarzschild schw(r_s);
    const double b_crit = std::sqrt(27.0)/2.0*r_s;

    // r_in = 1.2 r_s fica ABAIXO da esfera de fotons: sem isso os raios
    // quase-criticos enrolariam no vazio e o teto seria inatingivel.
    const PowerLawHalo halo(3.0, 1.0e8, 2.0, 1.2*r_s, 1.0e4*r_s);

    // =================================================================
    std::printf("\n  T39  INSENSIBILIDADE A NORMALIZACAO  (roda primeiro:\n");
    std::printf("       se R nao for invariante, nada a jusante significa nada)\n");
    {
        const std::string p = "/tmp/phasis_t39_ondulada.dat";
        escreve_tabela_ondulada(p);
        auto base = std::make_shared<TableCrossSection>(p);

        std::printf("    tabela ondulada: %zu pontos, E em [%.3g, %.3g]\n",
                    base->size(), base->E_min(), base->E_max());

        IntegratorOpts o;
        const double bs[3] = { b_crit*1.001, 3.0*r_s, 30.0*r_s };
        const double Es[3] = { 1.0e5, 1.0e7, 1.0e9 };

        double pior = 0.0;
        std::printf("    %-8s %-11s %-14s %-12s %-12s %s\n",
                    "k", "b/r_s", "E_inf", "R", "R(k=1)", "|dif| rel");
        for (int ib = 0; ib < 3; ++ib) {
            for (int ie = 0; ie < 3; ++ie) {
                const Rmed ref = medir_R(*base, schw, halo, bs[ib], Es[ie], o);
                for (double k : {0.5, 2.0, 10.0}) {
                    const ScaledCrossSection sc(base, k);
                    const Rmed m = medir_R(sc, schw, halo, bs[ib], Es[ie], o);
                    const double d = std::fabs(m.R - ref.R)/std::fabs(ref.R);
                    pior = std::max(pior, d);
                    if (ib == 0)
                        std::printf("    %-8.1f %-11.4f %-14.3g %-12.9f %-12.9f %.2e\n",
                                    k, bs[ib]/r_s, Es[ie], m.R, ref.R, d);
                }
            }
        }
        std::printf("    pior de 27 combinacoes (3 k x 3 b x 3 E): %.3e\n", pior);
        ++total;
        if (!(pior <= 1.0e-12)) { ++falhas;
            std::printf("    <-- FALHA: R depende da normalizacao\n"); }
        else std::printf("    OK: R invariante sob sigma -> k sigma em %.1e\n", pior);

        // Mesma coisa na tabela REAL do dipolo, se existir. A ondulada
        // e controlada; a real tem o ruido de quadratura de verdade.
        std::ifstream probe(kTabelaDipolo);
        if (probe) {
          try {
            auto dip = std::make_shared<TableCrossSection>(kTabelaDipolo);
            double pr = 0.0;
            for (double b : {b_crit*1.001, 5.0*r_s}) {
                for (double E : {1.0e7, 1.0e11}) {
                    if (E*std::sqrt(3.0) > dip->E_max()) continue;
                    const Rmed ref = medir_R(*dip, schw, halo, b, E, o);
                    for (double k : {0.5, 2.0, 10.0}) {
                        const ScaledCrossSection sc(dip, k);
                        const Rmed m = medir_R(sc, schw, halo, b, E, o);
                        pr = std::max(pr, std::fabs(m.R - ref.R)/std::fabs(ref.R));
                    }
                }
            }
            std::printf("    tabela real do dipolo: pior %.3e\n", pr);
            ++total;
            if (!(pr <= 1.0e-12)) { ++falhas;
                std::printf("    <-- FALHA: R depende da normalizacao (tabela real)\n"); }
            else std::printf("    OK: idem na tabela real\n");
          } catch (const std::exception& e) {
            // Uma tabela truncada ou malformada tem de virar FALHA
            // limpa, com a mensagem, nunca core dump: e o que acontece
            // se alguem ler o arquivo enquanto ele esta sendo escrito.
            ++total; ++falhas;
            std::printf("    <-- FALHA ao ler %s: %s\n", kTabelaDipolo, e.what());
          }
        } else {
            std::printf("    (tabela do dipolo ausente -- %s)\n", kTabelaDipolo);
        }
    }

    // =================================================================
    std::printf("\n  VALIDACAO CRUZADA de metadados entre duas tabelas\n");
    {
        // Este e o unico ponto onde uma comparacao errada passaria
        // despercebida. Cada tabela isolada e autoconsistente: uma de nu
        // e uma de nubar sao as duas secoes de choque validas, e a razao
        // entre elas continua sendo um numero. So a COMPARACAO e que
        // deixa de significar algo -- e nada dentro de cada arquivo
        // consegue detectar isso.
        const std::string a = "/tmp/phasis_meta_a.dat";
        const std::string b = "/tmp/phasis_meta_b.dat";
        escreve_com_meta(a, "nu", "isoscalar_nucleon", "CC");

        escreve_com_meta(b, "nu", "isoscalar_nucleon", "CC");
        chk_bool("nu x nu, mesmo alvo: passa",
                 lanca([&]{ assert_comparable(TableCrossSection(a),
                                              TableCrossSection(b)); }), false);

        escreve_com_meta(b, "nubar", "isoscalar_nucleon", "CC");
        chk_bool("nu x nubar: lanca", 
                 lanca([&]{ assert_comparable(TableCrossSection(a),
                                              TableCrossSection(b)); }), true);

        escreve_com_meta(b, "nu", "proton", "CC");
        chk_bool("isoescalar x proton: lanca",
                 lanca([&]{ assert_comparable(TableCrossSection(a),
                                              TableCrossSection(b)); }), true);

        escreve_com_meta(b, "nu", "isoscalar_nucleon", "NC");
        chk_bool("CC x NC: lanca",
                 lanca([&]{ assert_comparable(TableCrossSection(a),
                                              TableCrossSection(b)); }), true);

        escreve_com_meta(b, "nu", "isoscalar_nucleon", "CC", /*completa=*/false);
        chk_bool("faltando 3 das nove chaves: lanca",
                 lanca([&]{ assert_comparable(TableCrossSection(a),
                                              TableCrossSection(b)); }), true);

        // As duas tabelas de producao TEM de passar.
        std::ifstream p1(kTabelaDipolo), p2("data/sigma_nuN_CC_collinear.dat");
        if (p1 && p2) {
            chk_bool("dipolo x colinear de producao: passa",
                     lanca([&]{ assert_comparable(
                         TableCrossSection(kTabelaDipolo),
                         TableCrossSection("data/sigma_nuN_CC_collinear.dat")); }),
                     false);
        } else {
            std::printf("    (tabelas de producao ausentes -- pulando)\n");
        }
    }

    // =================================================================
    std::printf("\n  T36  CONSISTENCIA DO ESTIMADOR de d ln sigma / d ln E\n");
    {
        // Unico teste que valida o estimador. Trivial de proposito: com
        // lei de potencia a resposta e alpha, por construcao, em
        // qualquer E e qualquer passo.
        for (double alpha : {0.0, 0.2, 0.36, 1.0}) {
            PowerLawCrossSection xs(1.0e-31, 1.0e6, alpha);
            for (double E : {1.0e4, 1.0e6, 1.0e9, 1.0e12}) {
                const LogSlope s = xs.log_slope_detail(E);
                char nome[96];
                std::snprintf(nome, sizeof nome,
                              "alpha=%.2f em E=%.0e (h=%.3f, %d cortes)",
                              alpha, E, s.h, s.halvings);
                // alpha = 0 exige tolerancia ABSOLUTA: nao ha referencia
                // relativa a um zero.
                if (alpha == 0.0) {
                    ++total;
                    const bool ok = std::fabs(s.alpha) <= 1.0e-10;
                    if (!ok) ++falhas;
                    std::printf("    %-46s %-18.12g abs %8.2e            %s\n",
                                nome, s.alpha, std::fabs(s.alpha), ok ? "OK" : "<-- FALHA");
                } else {
                    chk(nome, s.alpha, alpha, 1.0e-10);
                }
            }
        }
        // O estimador tem de recuar para diferenca lateral na borda do
        // dominio de uma tabela, em vez de pedir sigma fora da faixa.
        const std::string p = "/tmp/phasis_t39_ondulada.dat";
        const TableCrossSection t(p);
        const LogSlope lo = t.log_slope_detail(t.E_min());
        const LogSlope hi = t.log_slope_detail(t.E_max());
        chk_bool("borda inferior usa diferenca lateral", lo.one_sided, true);
        chk_bool("borda superior usa diferenca lateral", hi.one_sided, true);
        chk_bool("meio da faixa usa diferenca centrada",
                 t.log_slope_detail(1.0e7).one_sided, false);
    }

    // =================================================================
    std::printf("\n  T37  TETO SATURADO:  max_b R -> 3^(alpha/2) - 1\n");
    {
        // A lei que torna isto afiado, derivada e depois confirmada
        // numericamente a 6 algarismos:
        //
        //   R = <g> - 1,  g = f^(-alpha/2),  peso w = n sigma dl
        //   teto - 1 - R = INT w (g_max - g) / INT w  =  K / tau_frozen
        //
        // O numerador CONVERGE quando b -> b_crit+, porque (g_max - g)
        // se anula exatamente onde o peso diverge. O denominador diverge
        // logaritmicamente. Logo o deficit cai como 1/tau_frozen, e dois
        // pontos bastam para extrapolar:
        //
        //   teto - 1 = (R1 tau1 - R2 tau2)/(tau1 - tau2)
        //
        // Nao e ajuste: e a solucao exata de duas equacoes numa reta.
        IntegratorOpts o; o.rel_tol = 1.0e-10;

        const double deltas[3] = { 1.0e-13, 1.0e-10, 1.0e-7 };
        for (double alpha : {0.2, 0.30, 0.36}) {
            const PowerLawCrossSection xs(1.0e-33, 1.0e6, alpha);
            const double teto = std::pow(3.0, 0.5*alpha);
            std::printf("    alpha = %.2f   3^(alpha/2) = %.9f\n", alpha, teto);

            double R[3], TAU[3], K[3];
            bool todos_ok = true, monotono = true, sob_o_teto = true;
            for (int i = 0; i < 3; ++i) {
                const double b = b_crit*(1.0 + deltas[i]);
                const Rmed m = medir_R(xs, schw, halo, b, 1.0e8, o);
                R[i] = m.R; TAU[i] = m.tau_frozen;
                K[i] = (teto - 1.0 - m.R)*m.tau_frozen;
                todos_ok = todos_ok && m.ok;
                // Cota EXATA por raio, sem limite nenhum: f(r) >= f(r_t)
                // ao longo de todo o caminho, logo <g> < g(r_t).
                const double cota = std::pow(schw.f(m.r_turn), -0.5*alpha) - 1.0;
                if (!(m.R < cota)) sob_o_teto = false;
                std::printf("      b/b_c-1=%-8.0e r_t/r_s=%-11.8f tau=%-11.4f "
                            "R=%-11.8f K=(teto-1-R)tau=%-10.5f\n",
                            deltas[i], m.r_turn/r_s, m.tau_frozen, m.R, K[i]);
            }
            for (int i = 0; i + 1 < 3; ++i) if (!(R[i] > R[i+1])) monotono = false;

            chk_bool("R cresce quando b -> b_crit+", monotono, true);
            chk_bool("R < f(r_t)^(-alpha/2) - 1 (cota exata por raio)",
                     sob_o_teto, true);
            chk_bool("quadratura satisfez a tolerancia", todos_ok, true);
            // K constante e o conteudo da lei; se ele nao for constante,
            // a extrapolacao seria um ajuste e nao uma identidade.
            chk("K constante entre b_c(1+1e-13) e b_c(1+1e-7)",
                K[0], K[2], 1.0e-5);

            const double ext = (R[0]*TAU[0] - R[1]*TAU[1])/(TAU[0] - TAU[1]);
            char nome[96];
            std::snprintf(nome, sizeof nome, "1+R extrapolado -> 3^(%.2f)", 0.5*alpha);
            chk(nome, 1.0 + ext, teto, 1.0e-8);
        }
        std::printf("    nota: os tres tetos sao 1.116123, 1.179148, 1.218658.\n");
        std::printf("          A especificacao dizia 1.128 e 1.176 para os dois\n");
        std::printf("          primeiros; 3^0.10 = 1.11612 e 3^0.15 = 1.17915.\n");
    }

    // =================================================================
    std::printf("\n  T38  MONOTONICIDADE DO TETO na tabela real (E >= 1e7 GeV)\n");
    {
        std::ifstream probe(kTabelaDipolo);
        if (!probe) {
            ++total; ++falhas;
            std::printf("    <-- FALHA: %s ausente. Gere com\n", kTabelaDipolo);
            std::printf("        cd dipole && make build/sigma_nuN && ./build/sigma_nuN ...\n");
        } else try {
            const TableCrossSection t(kTabelaDipolo);
            const double E0 = std::max(1.0e7, t.E_min()*1.05);
            const double E1 = t.E_max()*0.95;

            // (a) TENDENCIA, com linha de base de uma decada. E a
            //     afirmacao fisica, e ela e robusta ao ruido de no da
            //     tabela porque a linha de base e 27 vezes o
            //     espacamento entre nos.
            const int nd = 40;
            std::vector<double> Ed, ad;
            for (int i = 0; i < nd; ++i) {
                const double E = E0*std::pow(E1/E0, double(i)/(nd-1));
                const double lo = E/std::sqrt(10.0), hi = E*std::sqrt(10.0);
                if (lo < t.E_min() || hi > t.E_max()) continue;
                Ed.push_back(E);
                ad.push_back(std::log(t.sigma_tot(hi)/t.sigma_tot(lo))/std::log(10.0));
            }
            int sobe_d = 0; double pior_d = 0.0; double E_pior_d = 0.0;
            for (std::size_t i = 0; i + 1 < ad.size(); ++i) {
                const double d = ad[i+1] - ad[i];
                if (d > 0.0) { ++sobe_d; if (d > pior_d) { pior_d = d; E_pior_d = Ed[i+1]; } }
            }
            std::printf("    (a) tendencia, linha de base de uma decada:\n");
            std::printf("        alpha em 1e7=%.4f  1e9=%.4f  1e11=%.4f  1e13=%.4f\n",
                        ad.empty() ? 0.0 : ad.front(),
                        ad.size()>2 ? ad[ad.size()/3] : 0.0,
                        ad.size()>2 ? ad[2*ad.size()/3] : 0.0,
                        ad.empty() ? 0.0 : ad.back());
            std::printf("        teto de %.5f a %.5f;  segmentos que sobem: %d de %d\n",
                        std::pow(3.0, 0.5*(ad.empty()?0.0:ad.front())),
                        std::pow(3.0, 0.5*(ad.empty()?0.0:ad.back())),
                        sobe_d, int(ad.size()) - 1);
            if (sobe_d) std::printf("        maior subida %.2e em E = %.3g GeV\n",
                                    pior_d, E_pior_d);
            chk_bool("teto nao-crescente (linha de base de 1 decada)",
                     sobe_d == 0, true);

            // (b) PONTUAL, com o estimador tal como especificado. Aqui o
            //     ruido de no da tabela entra sem filtro, e o veredito
            //     tem de ser quantitativo: comparo cada subida com a
            //     incerteza que o proprio ruido da tabela injeta em
            //     alpha no passo h que o estimador escolheu.
            const int np = 90;
            std::vector<double> Ep, ap, hp;
            for (int i = 0; i < np; ++i) {
                const double E = E0*std::pow(E1/E0, double(i)/(np-1));
                const LogSlope s = t.log_slope_detail(E);
                Ep.push_back(E); ap.push_back(s.alpha); hp.push_back(s.h);
            }
            // ruido da tabela: residuo de ln sigma em relacao a um
            // ajuste quadratico local sobre 5 nos consecutivos.
            double soma = 0.0; int nres = 0;
            for (std::size_t i = 2; i + 2 < Ep.size(); ++i) {
                // segunda diferenca centrada de 5 pontos filtra a
                // tendencia suave e deixa o ruido
                const double v = ( t.sigma_tot(Ep[i-2]) > 0.0 ) ?
                    (std::log(t.sigma_tot(Ep[i-2])) - 4*std::log(t.sigma_tot(Ep[i-1]))
                     + 6*std::log(t.sigma_tot(Ep[i])) - 4*std::log(t.sigma_tot(Ep[i+1]))
                     + std::log(t.sigma_tot(Ep[i+2])))/std::sqrt(70.0) : 0.0;
                soma += v*v; ++nres;
            }
            const double eps = (nres > 0) ? std::sqrt(soma/nres) : 0.0;
            const double h_tip = hp.empty() ? 0.1 : hp[hp.size()/2];
            const double sigma_alpha = std::sqrt(2.0)*eps/(2.0*h_tip);

            int sobe_p = 0; double pior_p = 0.0, E_pior_p = 0.0;
            for (std::size_t i = 0; i + 1 < ap.size(); ++i) {
                const double d = ap[i+1] - ap[i];
                if (d > 0.0) { ++sobe_p; if (d > pior_p) { pior_p = d; E_pior_p = Ep[i+1]; } }
            }
            std::printf("    (b) pontual, estimador com h = %.4f em ln E:\n", h_tip);
            std::printf("        ruido de ln sigma na tabela: %.2e\n", eps);
            std::printf("        -> incerteza em alpha: %.4f ; queda esperada por passo: %.4f\n",
                        sigma_alpha,
                        ap.empty() ? 0.0 : (ap.front()-ap.back())/std::max<std::size_t>(1, ap.size()-1));
            std::printf("        subidas: %d de %d ; maior %.4f em E = %.3g GeV\n",
                        sobe_p, int(ap.size()) - 1, pior_p, E_pior_p);
            ++total;
            if (sobe_p == 0) {
                std::printf("    OK: teto pontual nao-crescente\n");
            } else if (pior_p <= 3.0*sigma_alpha) {
                std::printf("    OK: as subidas (max %.4f) cabem em 3 sigma do ruido\n"
                            "        da tabela (%.4f). O limite e da TABELA, nao do\n"
                            "        estimador: reduzir o ruido pede mais nos na\n"
                            "        quadratura de sigma no dipole.\n",
                            pior_p, 3.0*sigma_alpha);
            } else {
                ++falhas;
                std::printf("    <-- FALHA: subida de %.4f em E = %.3g GeV excede 3 sigma\n"
                            "        do ruido da tabela (%.4f). Isso nao e artefato de\n"
                            "        interpolacao: o problema esta em sigma(E).\n",
                            pior_p, E_pior_p, 3.0*sigma_alpha);
            }
        } catch (const std::exception& e) {
            ++total; ++falhas;
            std::printf("    <-- FALHA ao ler %s: %s\n", kTabelaDipolo, e.what());
        }
    }

    std::printf("\n=====================================================\n");
    std::printf("  %d verificacoes, %d falhas\n", total, falhas);
    std::printf("  RESULTADO: %s\n", falhas == 0 ? "OK" : "FALHA");
    std::printf("=====================================================\n\n");
    return falhas == 0 ? 0 : 1;
}
