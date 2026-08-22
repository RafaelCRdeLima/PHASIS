#ifndef PHASIS_CASCADE_HPP
#define PHASIS_CASCADE_HPP

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "phasis/density.hpp"
#include "phasis/integrate.hpp"
#include "phasis/metric.hpp"
#include "phasis/ray.hpp"
#include "phasis/cross_section.hpp"

namespace phasis {

// =====================================================================
// Regeneracao por corrente neutra.
//
// A variavel e E_inf, conservada ao longo da geodesica. O que torna isso
// tratavel:
//
//   Uma interacao NC em raio r leva E'_loc -> E_loc = (1-y) E'_loc. Como
//   o neutrino inicial e o final estao no MESMO r, o mesmo fator
//   sqrt(f(r)) aparece nos dois, e portanto
//
//       E_loc / E'_loc  =  E_inf / E'_inf  =  1 - y
//
//   ou seja, y e INVARIANTE SOB REDSHIFT. O mapeamento entre bins de
//   E_inf nao depende da posicao. So o COEFICIENTE varia ao longo do
//   raio, porque sigma e avaliada em E_loc -- a ESTRUTURA da matriz,
//   nunca.
//
// Consequencia arquitetural, e a razao de CascadeKernel existir: os
// intervalos de y que ligam bin j a bin i sao montados UMA VEZ, na
// construcao, e nunca mais. Remonta-los por passo seria O(N^2) por passo
// de EDO.
// =====================================================================

// ---------------------------------------------------------------------
struct DifferentialCrossSection {
    virtual ~DifferentialCrossSection() = default;

    virtual double dsigma_nc_dy(double E_GeV, double y) const = 0;
    virtual double sigma_nc(double E_GeV) const = 0;   // = integral em y
    virtual double sigma_cc(double E_GeV) const = 0;

    // Se true, dsigma/dy(E,y) = sigma_nc(E) * g(y) com g independente de
    // E. Nesse caso o kernel guarda G_ij = integral de g no intervalo, e
    // K_ij(E) = sigma_nc(E) * G_ij -- uma multiplicacao por passo, sem
    // nenhuma quadratura. Todos os kernels analiticos desta fase sao
    // assim. Um kernel de tabela com forma dependente de E nao e, e cai
    // no caminho tabelado.
    virtual bool shape_is_E_independent() const { return false; }

    // Pontos de quebra em y onde dsigma/dy nao e suave (bordas de
    // suporte, arestas de grade de tabela).
    //
    // O kernel PRECISA saber deles. A identidade de consistencia so
    // fecha se a soma dos pedacos for uma regra composta legitima sobre
    // [0,1]; uma descontinuidade DENTRO de um painel de Simpson quebra
    // isso, e o erro nao cancela entre intervalos adjacentes.
    virtual std::vector<double> y_breakpoints() const { return {}; }

    // Faixa em E onde dsigma/dy esta definida. A base nao restringe; a
    // tabela restringe a faixa tabelada. O kernel precisa disso para
    // saber ate onde pode tabelar G_ij(E) -- e para recusar, em vez de
    // extrapolar, quando o redshift empurra E_loc para fora.
    virtual double E_domain_min() const { return 0.0; }
    virtual double E_domain_max() const {
        return std::numeric_limits<double>::infinity();
    }

    virtual std::string name() const = 0;
};

// ---------------------------------------------------------------------
// dsigma_NC/dy = 0. Para T20: a cascata tem de reproduzir exp(-tau) das
// Fases 2-3 nos mesmos digitos.
class ZeroNC final : public DifferentialCrossSection {
public:
    ZeroNC(double sigma_cc0, double E0, double alpha)
        : s0_(sigma_cc0), E0_(E0), a_(alpha) {}
    double dsigma_nc_dy(double, double) const override { return 0.0; }
    double sigma_nc(double) const override { return 0.0; }
    double sigma_cc(double E) const override { return s0_*std::pow(E/E0_, a_); }
    bool shape_is_E_independent() const override { return true; }
    std::string name() const override { return "ZeroNC"; }
private:
    double s0_, E0_, a_;
};

// ---------------------------------------------------------------------
// dsigma_NC/dy uniforme em [0,1]. E o caso com SOLUCAO ANALITICA:
// g(y) = 1  =>  Z = integral_0^1 (1-y)^(gamma-1) dy = 1/gamma.
class ConstantY final : public DifferentialCrossSection {
public:
    ConstantY(double sigma_nc0, double sigma_cc0, double E0,
              double alpha_nc, double alpha_cc)
        : nc0_(sigma_nc0), cc0_(sigma_cc0), E0_(E0),
          anc_(alpha_nc), acc_(alpha_cc) {}

    double dsigma_nc_dy(double E, double y) const override {
        return (y >= 0.0 && y <= 1.0) ? sigma_nc(E) : 0.0;
    }
    double sigma_nc(double E) const override { return nc0_*std::pow(E/E0_, anc_); }
    double sigma_cc(double E) const override { return cc0_*std::pow(E/E0_, acc_); }
    bool shape_is_E_independent() const override { return true; }
    std::string name() const override { return "ConstantY"; }

    // Z analitico exato para fluxo E^-gamma.
    static double Z_analytic(double gamma) { return 1.0/gamma; }
private:
    double nc0_, cc0_, E0_, anc_, acc_;
};

// ---------------------------------------------------------------------
// dsigma_NC/dy proporcional a y^(-beta) em [y_min, 1], normalizada a
// sigma_nc(E). Kernel PICADO -- existe para T24 morder: grade log
// uniforme resolve mal o pico em y_min, e a ordem de convergencia deve
// degradar em relacao ao caso suave.
class PowerLawY final : public DifferentialCrossSection {
public:
    PowerLawY(double sigma_nc0, double sigma_cc0, double E0,
              double alpha_nc, double alpha_cc,
              double beta, double y_min);

    double dsigma_nc_dy(double E, double y) const override;
    double sigma_nc(double E) const override { return nc0_*std::pow(E/E0_, anc_); }
    double sigma_cc(double E) const override { return cc0_*std::pow(E/E0_, acc_); }
    bool shape_is_E_independent() const override { return true; }
    std::vector<double> y_breakpoints() const override { return { y_min_ }; }
    std::string name() const override { return "PowerLawY"; }

    // Z analitico: (1/norm) * integral_{y_min}^1 y^-beta (1-y)^(gamma-1) dy
    double Z_analytic(double gamma, int n_quad = 2000000) const;

    double beta()  const { return beta_; }
    double y_min() const { return y_min_; }
private:
    double nc0_, cc0_, E0_, anc_, acc_, beta_, y_min_, norm_;
};

// ---------------------------------------------------------------------
// Grade logaritmica em E_inf. Bins de largura igual em ln E.
class EnergyGrid {
public:
    EnergyGrid(double E_min, double E_max, int n_bins);

    int    n()        const { return n_; }
    double edge(int i) const;          // borda inferior do bin i; edge(n) = E_max
    double center(int i) const;        // media geometrica das bordas
    double ratio()    const { return ratio_; }
    double E_min()    const { return E_min_; }
    double E_max()    const { return E_max_; }

    // Integral de uma lei de potencia A E^-gamma sobre o bin i.
    double powerlaw_bin_integral(int i, double gamma) const;

private:
    double E_min_, E_max_, ratio_, ln_min_, dln_;
    int n_;
};

// ---------------------------------------------------------------------
// Kernel de acoplamento entre bins.
//
// Um neutrino no bin j (energia representativa E_j) que sofre NC com
// inelasticidade y vai para E = (1-y) E_j. O bin i recebe o intervalo
//
//     y em [ 1 - edge(i+1)/E_j ,  1 - edge(i)/E_j ]
//
// intersectado com [0,1]. Esses intervalos, para i = j, j-1, ..., 0,
// mais o vazamento (y acima de 1 - E_min/E_j), PARTICIONAM [0,1]
// exatamente. Dai a identidade de consistencia
//
//     soma_i K_ij  +  leakage_j  ==  sigma_nc(E_j)
//
// que e exigida na construcao. Se ela nao fecha, a cascata nao conserva
// numero e nada a jusante vale.
//
// O vazamento e RASTREADO, nunca descartado nem absorvido no bin de
// baixo -- absorver empilharia fluxo artificialmente na borda.
// FORMA DEPENDENTE DE E -- o caminho que o codigo chamava de Fase 4b.
//
// Uma tabela de dipolo real NAO tem dsigma/dy = sigma_nc(E) g(y) com g
// fixo: o propagador e a faixa de x acessivel mudam com a energia, e
// com eles a forma em y. Entao G_ij tem de virar G_ij(E).
//
// A rota escolhida: tabelar G_ij em nos LOG-E e interpolar LINEARMENTE
// em ln E. Isso nao e so conveniencia -- e o que preserva a fisica:
//
//     A identidade de consistencia  soma_i G_ij + G_leak_j = 1
//     e LINEAR nos G. Interpolacao linear de um conjunto que soma 1 em
//     cada no ainda soma 1 entre os nos.
//
// Ou seja: a cascata conserva numero EXATAMENTE em qualquer E, nao so
// nos nos da tabela. Com interpolacao cubica isso deixaria de valer, e
// o erro de conservacao apareceria como fluxo criado ou destruido do
// nada ao longo do raio -- que e justamente o tipo de erro que nenhuma
// checagem a jusante pegaria.
//
// Os intervalos em y que ligam o bin j ao bin i NAO dependem de E_loc:
// y e invariante sob redshift, e as bordas dos bins estao em E_inf. So
// o ARGUMENTO de dsigma/dy varia. Por isso da para tabelar em E de uma
// vez e nunca mais reintegrar.
class CascadeKernel {
public:
    CascadeKernel(const EnergyGrid& grid,
                  const DifferentialCrossSection& xsec,
                  double consistency_tol = 1.0e-12,
                  int n_y_quad = 4096,
                  // Nos em ln E, usados so quando a forma depende de E.
                  int n_E_nodes = 24,
                  // Faixa em E a tabelar. 0 = automatica:
                  // [E_min/2, E_max*2] da grade, cortada pelo dominio da
                  // secao de choque. O fator 2 cobre com folga o boost
                  // maximo sqrt(3) do redshift em Schwarzschild.
                  double E_lo = 0.0,
                  double E_hi = 0.0);

    // Forma normalizada: G_ij = integral de g(y) no intervalo (i,j).
    // Vale K_ij(E) = sigma_nc(E) * G_ij(E).
    //
    // A versao sem E so e legitima quando a forma independe de E, e
    // LANCA caso contrario -- silenciosamente devolver o primeiro no
    // seria dar um numero plausivel e errado.
    double G(int i, int j) const;
    double G_leak(int j)   const;

    double G(int i, int j, double E) const;
    double G_leak(int j, double E)   const;

    bool shape_depends_on_E() const { return nE_ > 1; }
    double kernel_E_min() const { return E_nodes_.empty() ? 0.0 : std::exp(lnE_.front()); }
    double kernel_E_max() const { return E_nodes_.empty() ? 0.0 : std::exp(lnE_.back()); }
    int n_E_nodes() const { return nE_; }

    // Residuo da identidade discreta, por bin j.
    // Residuo da identidade discreta, por bin j (pior sobre os nos de E).
    double consistency_residual(int j) const;

    // Quanto a soma dos pedacos, calculada com o integrador adaptativo,
    // difere de sigma_nc(E) da propria secao de choque. Num kernel
    // analitico e zero; num kernel de tabela mede a RESOLUCAO da grade
    // em y -- a tabela integra por trapezio sobre os nos, o kernel
    // integra o interpolante de forma adaptativa, e a diferenca e o
    // erro do trapezio.
    double worst_norm_mismatch() const;

    // true quando sigma_nc == 0: nao ha regeneracao, e a identidade
    // normalizada nao se aplica (a real, soma K + leak = 0, e trivial).
    bool no_regeneration() const { return no_regeneration_; }
    double worst_residual() const;

    const EnergyGrid& grid() const { return grid_; }
    const DifferentialCrossSection& xsec() const { return xsec_; }

    // Fator Z da solucao analitica, na forma DISCRETA:
    //     Z_disc(gamma) = soma_d G_d * rho^(-d(gamma-1))
    // Ver cascade.cpp para a derivacao e para o que ele significa.
    double Z_discrete(double gamma) const;

private:
    std::size_t idx(int e, int i, int j) const {
        return (static_cast<std::size_t>(e)*static_cast<std::size_t>(n_)
                + static_cast<std::size_t>(j))*static_cast<std::size_t>(n_)
               + static_cast<std::size_t>(i);
    }
    // peso e nos vizinhos em ln E
    void locate(double E, int& e0, int& e1, double& t) const;

    const EnergyGrid& grid_;
    const DifferentialCrossSection& xsec_;
    int n_;
    int nE_ = 1;
    std::vector<double> E_nodes_;
    std::vector<double> lnE_;
    std::vector<double> G_;        // nE_ x n_ x n_, so i <= j e nao nulo
    std::vector<double> G_leak_;   // nE_ x n_
    std::vector<double> residual_; // nE_ x n_
    std::vector<double> norm_mismatch_;
    bool no_regeneration_ = false;
};

// ---------------------------------------------------------------------
// dsigma_NC/dy lida de tabela, com METADADOS OBRIGATORIOS.
//
// Por que obrigatorios: nada neste projeto forca a explicitar a
// convencao, porque quem escreve e quem le sao o mesmo autor. E a
// convencao e onde mora o fator 2 -- ja mordeu uma vez, no acoplamento
// NC do dipole. O leitor EXIGE cada chave e valida `convention_y` e
// `current` contra o que o chamador declara esperar. Ausencia de
// qualquer uma e excecao, nunca default.
//
// Grade log(E) x log(y), interpolacao bilinear em log-log. Sem
// extrapolacao. Uma lei de potencia em qualquer dos eixos e reta neste
// plano, logo e reproduzida EXATAMENTE -- e o que o round-trip verifica.
class TableDifferentialCrossSection final : public DifferentialCrossSection {
public:
    struct Expect {
        std::string convention_y = "(E_in - E_out)/E_in";
        std::string current      = "NC";
    };

    TableDifferentialCrossSection(const std::string& path,
                                  const Expect& esperado,
                                  std::shared_ptr<const CrossSection> sigma_cc);

    double dsigma_nc_dy(double E_GeV, double y) const override;
    double sigma_nc(double E_GeV) const override;
    double sigma_cc(double E_GeV) const override;
    std::string name() const override { return "TableDifferential[" + path_ + "]"; }

    double E_domain_min() const override { return std::exp(lnE_.front()); }
    double E_domain_max() const override { return std::exp(lnE_.back()); }

    // Bordas da grade em y da tabela. Sao pontos de quebra reais: o
    // corte cinematico y < Q2min/(s x_max) poe um degrau em y, e uma
    // descontinuidade DENTRO de um painel de Simpson nao cancela entre
    // intervalos vizinhos.
    std::vector<double> y_breakpoints() const override;

    double y_min_table() const { return std::exp(lny_.front()); }

    const std::string& meta(const std::string& chave) const;

private:
    std::string path_;
    std::vector<double> lnE_, lny_;
    std::vector<double> lnD_;      // ln(dsigma/dy), nE x nY
    std::vector<double> snc_;      // sigma_nc por no de E
    std::shared_ptr<const CrossSection> cc_;
    std::vector<std::pair<std::string,std::string>> meta_;
};

// Escreve uma tabela no formato acima a partir de um kernel analitico.
// Existe para o round-trip: valida o LEITOR sem depender do dipole.
void write_dsigma_table(const std::string& path,
                        const DifferentialCrossSection& xsec,
                        double E_min, double E_max, int nE,
                        double y_min, int nY,
                        const std::string& gerado_por);

// Z_disc calculado SO para o bin do topo, sem montar a matriz N x N.
//
// T24 precisa varrer bins/decada por duas decadas, e a construcao
// completa do kernel e O(N^2) integrais -- inviavel em N ~ 4000. Esta
// rota e O(N).
double z_discrete_row(const EnergyGrid& grid,
                      const DifferentialCrossSection& xsec,
                      double gamma);

// =====================================================================
// Transporte do espectro ao longo de uma geodesica.
//
//   dphi_i/dl = -n(l) sigma_tot(E_loc_i) phi_i
//               + n(l) soma_{j>=i} K_ij(E_loc_j) phi_j
//   dleak/dl  =  n(l) soma_j  Kleak_j(E_loc_j) phi_j
//
// com K_ij(E) = sigma_nc(E) G_ij  --  a ESTRUTURA G_ij nao depende de l,
// so o coeficiente sigma_nc(E_loc). E o que torna o passo O(N^2) em
// multiplicacoes, sem nenhuma quadratura por passo.
//
// O vazamento e acumulado numa componente PROPRIA do estado, nunca
// descartado nem empilhado no bin de baixo.
// =====================================================================
struct CascadeResult {
    std::vector<double> phi;        // espectro na saida
    double leakage_total = 0.0;
    double column_g_cm2  = 0.0;
    bool   tolerance_met = true;
    bool   captured      = false;
    long   n_evals       = 0;
};

CascadeResult transport_cascade(const std::vector<double>& phi0,
                                const CascadeKernel& kernel,
                                const Ray& ray,
                                const Metric& metric,
                                const DensityProfile& profile,
                                const IntegratorOpts& opts,
                                // T25: congela a matriz em E_inf, ignorando o
                                // redshift no argumento de sigma. Mantem a
                                // GEOMETRIA identica -- comparar contra
                                // Minkowski mudaria o caminho tambem, e a
                                // diferenca medida nao seria so do redshift.
                                bool freeze_redshift = false);

} // namespace phasis

#endif
