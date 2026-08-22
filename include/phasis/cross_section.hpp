#ifndef PHASIS_CROSS_SECTION_HPP
#define PHASIS_CROSS_SECTION_HPP

#include <cmath>
#include <iosfwd>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace phasis {

// =====================================================================
// Resultado do estimador de d ln sigma / d ln E, com o diagnostico que
// permite julga-lo.
//
// ESCOLHA DO ESTIMADOR -- diferenca centrada, nao derivada analitica.
//
// A alternativa analitica so seria legitima se a interpolacao fosse
// C^1. Nao e: TableCrossSection interpola LINEARMENTE em (ln E, ln
// sigma), o que e C^0. A derivada analitica dessa spline e uma funcao
// ESCADA, com salto em cada no, e no proprio no nem esta definida.
// Entao a rota analitica esta fechada por construcao, e a centrada e a
// unica honesta.
//
// O passo e adaptativo por REDUCAO PELA METADE, e converge por um
// motivo especifico da forma da tabela: enquanto a janela [E e^-h,
// E e^+h] contem um no, D(h) e uma media ponderada de dois segmentos e
// muda ao encolher; assim que a janela cabe DENTRO de um segmento,
// D(h) para de mudar porque ali ln sigma e exatamente linear em ln E.
// O criterio de parada e portanto "encontrei o segmento", nao "reduzi o
// erro de truncamento".
//
// Piso em h: sem ele o cancelamento catastrofico domina. Os valores de
// ln sigma valem ~ -70 (sigma ~ 1e-31 cm^2), com ulp ~ 1.4e-14, e a
// diferenca vale 2 h alpha. Com h = 3e-3 e alpha = 0.2 isso da erro
// relativo ~ 2e-11 -- ja no limite do que T36 exige. Por isso o
// estimador comeca GRANDE (h = 0.1) e so encolhe quando precisa: para
// uma lei de potencia exata o truncamento e nulo em qualquer h, e o h
// grande e estritamente melhor.
struct LogSlope {
    double alpha     = 0.0;
    double h         = 0.0;      // passo em ln E efetivamente usado
    int    halvings  = 0;
    bool   one_sided = false;    // E colado numa borda do dominio
    bool   converged = true;     // false: nao estabilizou ate o piso de h
};

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

    // Faixa em que sigma esta DEFINIDA. A base nao restringe; a tabela
    // restringe a faixa tabelada. log_slope precisa disso: sem saber a
    // borda, a diferenca centrada pediria sigma fora da faixa e a
    // tabela lancaria -- que e o comportamento correto dela.
    virtual double E_domain_min() const { return 0.0; }
    virtual double E_domain_max() const {
        return std::numeric_limits<double>::infinity();
    }

    // alpha_eff(E) = d ln sigma / d ln E.
    //
    // NAO e virtual de proposito. Se PowerLawCrossSection devolvesse
    // alpha analiticamente, T36 -- o unico teste que valida o estimador
    // -- estaria validando um return de campo. O estimador tem de ser o
    // MESMO nos dois casos.
    LogSlope log_slope_detail(double E_GeV) const;
    double   log_slope(double E_GeV) const { return log_slope_detail(E_GeV).alpha; }

    // Teto do acoplamento GR x saturacao para esta secao de choque em E:
    //
    //     3^(alpha_eff(E)/2)
    //
    // Vem de E_loc/E_inf = 1/sqrt(f(r_t)) <= sqrt(3) em Schwarzschild,
    // porque todo raio nao capturado tem r_t >= 1.5 r_s. Com sigma ~
    // E^alpha a razao de secoes de choque e no maximo 3^(alpha/2).
    double gr_saturation_ceiling(double E_GeV) const {
        return std::pow(3.0, 0.5*log_slope(E_GeV));
    }
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

    double E_domain_min() const override { return E_.front(); }
    double E_domain_max() const override { return E_.back(); }

    const std::string& path() const { return path_; }

    // Cabecalho '#' preservado, para poder ecoar a proveniencia da
    // tabela do dipole na saida do PHASIS.
    const std::vector<std::string>& header() const { return header_; }

    // --- metadados -------------------------------------------------
    // Linhas '# chave = valor' viram pares. O formato antigo do dipole
    // ('# chave valor', sem '='), continua sendo aceito como cabecalho
    // livre e simplesmente nao vira metadado -- por isso o construtor
    // NAO exige nada. Quem exige e assert_comparable, no unico ponto
    // onde a ausencia faria diferenca.
    bool has_meta(const std::string& chave) const;
    const std::string& meta(const std::string& chave) const;   // lanca se faltar
    const std::vector<std::pair<std::string,std::string>>& metadata() const {
        return meta_;
    }

    // Lanca se qualquer uma das nove chaves obrigatorias faltar.
    void require_full_metadata() const;

private:
    std::string path_;
    std::vector<double> E_;      // crescente
    std::vector<double> sigma_;
    std::vector<double> lnE_;
    std::vector<double> lnSigma_;
    std::vector<std::string> header_;
    std::vector<std::pair<std::string,std::string>> meta_;
};

// As nove chaves que toda tabela de secao de choque deve carregar.
// Mesma lista que cascade.cpp exige da tabela de dsigma/dy.
const std::vector<std::string>& required_metadata_keys();

// ---------------------------------------------------------------------
// Validacao CRUZADA entre duas tabelas carregadas ao mesmo tempo.
//
// Este e o unico ponto do projeto onde uma comparacao errada passaria
// despercebida. Cada tabela isolada e autoconsistente: se uma for de
// nubar e a outra de nu, ou uma de alvo isoescalar e a outra de proton,
// as duas continuam sendo secoes de choque perfeitamente validas, e a
// razao entre elas continua sendo um numero. So a comparacao e que
// deixa de significar alguma coisa -- e nada dentro de cada tabela
// consegue detectar isso.
//
// Exige as nove chaves nas duas e IGUALDADE de target, projectile,
// current, units_sigma e units_E. Nao exige igualdade de dipole_model
// nem de generated_by: e exatamente ali que as duas TEM de diferir.
void assert_comparable(const TableCrossSection& a, const TableCrossSection& b);

// Ecoa os metadados de uma tabela num cabecalho de CSV, cada linha
// prefixada por `tag`. Existe para que um arquivo de saida denuncie
// sozinho, daqui a seis meses, se alguem comparou maca com laranja.
void write_metadata_header(std::ostream& os,
                           const std::string& tag,
                           const TableCrossSection& t);

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

    // O dominio e o da base. Sem isto, log_slope de uma tabela
    // reescalada usaria o dominio infinito da base abstrata, pediria
    // sigma fora da faixa e a tabela lancaria -- e T39, que compara a
    // reescalada com a original, quebraria por um motivo que nao tem
    // nada a ver com o que ele mede.
    double E_domain_min() const override { return base_->E_domain_min(); }
    double E_domain_max() const override { return base_->E_domain_max(); }

    double k() const { return k_; }
    const CrossSection& base() const { return *base_; }

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

    // Interseccao dos dominios: uma soma so esta definida onde TODAS as
    // parcelas estao.
    double E_domain_min() const override;
    double E_domain_max() const override;

private:
    std::vector<std::shared_ptr<const CrossSection>> parts_;
};

} // namespace phasis

#endif
