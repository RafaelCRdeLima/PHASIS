#ifndef PHASIS_CASCADE_HPP
#define PHASIS_CASCADE_HPP

#include <cmath>
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
class CascadeKernel {
public:
    CascadeKernel(const EnergyGrid& grid,
                  const DifferentialCrossSection& xsec,
                  double consistency_tol = 1.0e-12,
                  int n_y_quad = 4096);

    // Forma normalizada: G_ij = integral de g(y) no intervalo (i,j).
    // Vale K_ij(E) = sigma_nc(E) * G_ij quando a forma independe de E.
    double G(int i, int j) const { return G_[idx(i, j)]; }
    double G_leak(int j)   const { return G_leak_[static_cast<std::size_t>(j)]; }

    // Residuo da identidade discreta, por bin j.
    double consistency_residual(int j) const { return residual_[static_cast<std::size_t>(j)]; }

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
    std::size_t idx(int i, int j) const {
        return static_cast<std::size_t>(j)*static_cast<std::size_t>(n_) + static_cast<std::size_t>(i);
    }

    const EnergyGrid& grid_;
    const DifferentialCrossSection& xsec_;
    int n_;
    std::vector<double> G_;        // n_ x n_, so i <= j e nao nulo
    std::vector<double> G_leak_;
    std::vector<double> residual_;
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
