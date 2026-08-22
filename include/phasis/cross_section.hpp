#ifndef PHASIS_CROSS_SECTION_HPP
#define PHASIS_CROSS_SECTION_HPP

#include <memory>
#include <string>
#include <vector>

namespace phasis {

// =====================================================================
// Secao de choque total neutrino-nucleon, em cm^2, funcao da energia
// em GeV.
//
// A energia pedida e sempre a LOCAL (medida pelo observador estatico no
// ponto), nunca a conservada no infinito. Na Fase 2 as duas diferem pelo
// redshift, e sigma passa a variar ao longo do raio.
// =====================================================================
struct CrossSection {
    virtual ~CrossSection() = default;
    virtual double sigma_tot(double E_GeV) const = 0;
    virtual std::string name() const = 0;
};

// ---------------------------------------------------------------------
// sigma = sigma0 * (E/E0)^alpha.  Existe para os testes: com ela varias
// integrais tem forma fechada.
class PowerLawCrossSection final : public CrossSection {
public:
    PowerLawCrossSection(double sigma0, double E0, double alpha)
        : sigma0_(sigma0), E0_(E0), alpha_(alpha) {}

    double sigma_tot(double E_GeV) const override;
    std::string name() const override { return "PowerLawCrossSection"; }

private:
    double sigma0_, E0_, alpha_;
};

// ---------------------------------------------------------------------
// Tabela lida de arquivo, interpolada LINEARMENTE EM LOG-LOG.
//
// Formatos aceitos (detectados por linha):
//
//   a) CSV de duas colunas:      E_GeV,sigma_cm2
//   b) Saida do dipole (3 col):  Enu_GeV  sigma_GeV_minus2  sigma_cm2
//
// O formato (b) e o que PHASIS/dipole produz; as colunas 1 e 3 sao
// usadas. Linhas comecando com '#' e linhas vazias sao ignoradas, o que
// preserva o cabecalho de proveniencia das tabelas do dipole.
//
// NAO EXTRAPOLA. Fora da faixa tabelada lanca std::out_of_range. Isso e
// deliberado: extrapolar uma secao de choque de saturacao para fora do
// intervalo em que ela foi calculada e exatamente o tipo de coisa que
// produz resultado errado sem aviso.
class TableCrossSection final : public CrossSection {
public:
    explicit TableCrossSection(const std::string& path);

    double sigma_tot(double E_GeV) const override;
    std::string name() const override { return "TableCrossSection[" + path_ + "]"; }

    double E_min() const { return E_[0]; }
    double E_max() const { return E_.back(); }
    std::size_t size() const { return E_.size(); }

    // Cabecalho '#' preservado, para poder ecoar a proveniencia da
    // tabela do dipole na saida do PHASIS.
    const std::vector<std::string>& header() const { return header_; }

private:
    std::string path_;
    std::vector<double> E_;      // crescente
    std::vector<double> sigma_;
    std::vector<double> lnE_;
    std::vector<double> lnSigma_;
    std::vector<std::string> header_;
};

// ---------------------------------------------------------------------
// sigma = k * base.
//
// Existe para uma coisa concreta: as tabelas produzidas por
// PHASIS/dipole sao de CORRENTE CARREGADA APENAS. Somar a corrente
// neutra exige ou uma segunda tabela, ou uma razao sigma_NC/sigma_CC
// declarada explicitamente.
class ScaledCrossSection final : public CrossSection {
public:
    ScaledCrossSection(std::shared_ptr<const CrossSection> base, double k)
        : base_(std::move(base)), k_(k) {}

    double sigma_tot(double E_GeV) const override {
        return k_*base_->sigma_tot(E_GeV);
    }
    std::string name() const override;

private:
    std::shared_ptr<const CrossSection> base_;
    double k_;
};

// ---------------------------------------------------------------------
// sigma = soma das partes. Para combinar CC + NC vindos de tabelas
// separadas.
class SumCrossSection final : public CrossSection {
public:
    void add(std::shared_ptr<const CrossSection> p) { parts_.push_back(std::move(p)); }

    double sigma_tot(double E_GeV) const override;
    std::string name() const override;

private:
    std::vector<std::shared_ptr<const CrossSection>> parts_;
};

} // namespace phasis

#endif
