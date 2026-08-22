#include "phasis/cross_section.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace phasis {

// =====================================================================
// alpha_eff(E) = d ln sigma / d ln E, por diferenca centrada em ln E
// com passo adaptativo.
//
// Ver o comentario de LogSlope em cross_section.hpp para POR QUE a rota
// analitica esta fechada (a interpolacao da tabela e C^0, nao C^1).
//
// Detalhe que importa: o denominador NAO e 2h, e sim
// ln(E_mais) - ln(E_menos), calculado nos MESMOS valores que sigma
// recebeu. E^h e depois ln nao devolvem h exatamente, e usar 2h
// injetaria um erro relativo ~1e-16/h -- em h = 3e-3 isso e 3e-14, na
// mesma ordem do que T36 mede.
namespace {

constexpr double kH0     = 0.1;      // passo inicial em ln E
constexpr double kHMin   = 3.0e-3;   // piso: abaixo disso o cancelamento manda
constexpr int    kMaxHalv = 8;
constexpr double kStable = 1.0e-11;  // "D(h) parou de mudar"

double slope_at(const CrossSection& xs, double E, double h, bool one_sided_up,
                bool one_sided_down)
{
    double Ep, Em;
    if (one_sided_up)        { Ep = E*std::exp(h); Em = E; }
    else if (one_sided_down) { Ep = E;             Em = E*std::exp(-h); }
    else                     { Ep = E*std::exp(h); Em = E*std::exp(-h); }

    const double den = std::log(Ep) - std::log(Em);
    if (!(den > 0.0)) throw std::domain_error("log_slope: passo degenerado");

    const double sp = xs.sigma_tot(Ep);
    const double sm = xs.sigma_tot(Em);
    if (!(sp > 0.0) || !(sm > 0.0)) {
        throw std::domain_error("log_slope: sigma <= 0, ln nao definido");
    }
    return (std::log(sp) - std::log(sm))/den;
}

} // namespace

LogSlope CrossSection::log_slope_detail(double E_GeV) const
{
    if (!(E_GeV > 0.0)) throw std::domain_error("log_slope: E <= 0");

    const double lo = E_domain_min();
    const double hi = E_domain_max();

    if (E_GeV < lo || E_GeV > hi) {
        std::ostringstream m;
        m << "log_slope: E = " << E_GeV << " fora do dominio ["
          << lo << ", " << hi << "] de " << name();
        throw std::out_of_range(m.str());
    }

    // Espaco disponivel de cada lado, em ln E. O fator 0.999 evita que
    // o ponto de avaliacao caia EXATAMENTE na borda, onde o
    // arredondamento de exp/log pode empurra-lo um ulp para fora e a
    // tabela -- corretamente -- lancar.
    const double room_lo = (lo > 0.0)
        ? 0.999*(std::log(E_GeV) - std::log(lo))
        : std::numeric_limits<double>::infinity();
    const double room_hi = std::isfinite(hi)
        ? 0.999*(std::log(hi) - std::log(E_GeV))
        : std::numeric_limits<double>::infinity();

    LogSlope out;

    bool up = false, down = false;
    double h = kH0;

    if (room_lo >= kHMin && room_hi >= kHMin) {
        h = std::min(kH0, std::min(room_lo, room_hi));
    } else if (room_hi >= kHMin) {
        up = true;  out.one_sided = true;
        h = std::min(kH0, room_hi);
    } else if (room_lo >= kHMin) {
        down = true; out.one_sided = true;
        h = std::min(kH0, room_lo);
    } else {
        std::ostringstream m;
        m << "log_slope: dominio estreito demais em torno de E = " << E_GeV
          << " para formar uma diferenca (precisa de " << kHMin << " em ln E)";
        throw std::domain_error(m.str());
    }

    double D = slope_at(*this, E_GeV, h, up, down);
    out.converged = false;

    for (int k = 0; k < kMaxHalv; ++k) {
        const double h2 = 0.5*h;
        if (h2 < kHMin) break;
        const double D2 = slope_at(*this, E_GeV, h2, up, down);
        ++out.halvings;
        if (std::fabs(D2 - D) <= kStable*std::max(1.0, std::fabs(D2))) {
            // Estabilizou. Devolve o h MAIOR: mesmo valor, menos
            // cancelamento.
            out.converged = true;
            break;
        }
        h = h2;
        D = D2;
    }

    out.alpha = D;
    out.h     = h;
    return out;
}

// ---------------------------------------------------------------------

double PowerLawCrossSection::sigma_tot(double E_GeV) const
{
    if (!(E_GeV > 0.0)) {
        throw std::domain_error("PowerLawCrossSection: E <= 0");
    }
    return sigma0_*std::pow(E_GeV/E0_, alpha_);
}

// ---------------------------------------------------------------------

TableCrossSection::TableCrossSection(const std::string& path) : path_(path)
{
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("TableCrossSection: nao consegui abrir " + path);
    }

    std::string line;
    while (std::getline(in, line)) {

        // preserva o cabecalho: e onde o dipole grava a proveniencia
        if (!line.empty() && line[0] == '#') {
            header_.push_back(line);
            const auto eq = line.find('=');
            if (eq != std::string::npos) {
                auto trim = [](std::string t) {
                    const auto a = t.find_first_not_of(" \t#\r\n");
                    if (a == std::string::npos) return std::string{};
                    const auto b = t.find_last_not_of(" \t\r\n");
                    return t.substr(a, b - a + 1);
                };
                const std::string chave = trim(line.substr(0, eq));
                // Chave tem de ser UM identificador. Sem esta regra,
                // qualquer linha de prosa com um '=' no meio vira
                // metadado: a linha
                //   '# largeXFactor (1-x)^7  [contagem, n_s=4]'
                //   virava a chave 'largeXFactor (1-x)^7 [contagem, n_s'
                // com valor '4', e ia parar no cabecalho dos CSV de
                // saida como se fosse proveniencia de verdade.
                const bool identificador =
                    !chave.empty() &&
                    chave.find_first_of(" \t") == std::string::npos;
                if (identificador) {
                    meta_.emplace_back(chave, trim(line.substr(eq + 1)));
                }
            }
            continue;
        }
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

        double E = 0.0, s = 0.0;
        bool ok = false;

        if (line.find(',') != std::string::npos) {
            // formato (a): CSV  E_GeV,sigma_cm2
            std::string tmp = line;
            std::replace(tmp.begin(), tmp.end(), ',', ' ');
            std::istringstream ss(tmp);
            ok = static_cast<bool>(ss >> E >> s);
        } else {
            // formato (b): saida do dipole
            //   Enu_GeV  sigma_GeV_minus2  sigma_cm2
            // ou duas colunas separadas por espaco.
            std::istringstream ss(line);
            double c1 = 0.0, c2 = 0.0, c3 = 0.0;
            if (ss >> c1 >> c2) {
                E = c1;
                if (ss >> c3) s = c3;   // 3 colunas: sigma em cm^2 e a terceira
                else          s = c2;   // 2 colunas
                ok = true;
            }
        }

        if (!ok) continue;
        if (!(E > 0.0) || !(s > 0.0)) continue;

        E_.push_back(E);
        sigma_.push_back(s);
    }

    if (E_.size() < 2) {
        throw std::runtime_error("TableCrossSection: menos de 2 pontos validos em " + path);
    }
    for (std::size_t i = 1; i < E_.size(); ++i) {
        if (!(E_[i] > E_[i-1])) {
            throw std::runtime_error(
                "TableCrossSection: grade de energia nao estritamente crescente em " + path);
        }
    }

    lnE_.reserve(E_.size());
    lnSigma_.reserve(E_.size());
    for (std::size_t i = 0; i < E_.size(); ++i) {
        lnE_.push_back(std::log(E_[i]));
        lnSigma_.push_back(std::log(sigma_[i]));
    }
}

double TableCrossSection::sigma_tot(double E_GeV) const
{
    if (!(E_GeV > 0.0)) {
        throw std::domain_error("TableCrossSection: E <= 0");
    }

    // Sem extrapolacao, por decisao de projeto. Estender uma secao de
    // choque de saturacao para fora da faixa em que foi calculada
    // produz resultado errado silenciosamente.
    if (E_GeV < E_.front() || E_GeV > E_.back()) {
        std::ostringstream m;
        m << "TableCrossSection: E = " << E_GeV << " GeV fora da faixa tabelada ["
          << E_.front() << ", " << E_.back() << "] em " << path_;
        throw std::out_of_range(m.str());
    }

    const double x = std::log(E_GeV);

    const auto it = std::upper_bound(lnE_.begin(), lnE_.end(), x);
    std::size_t j = static_cast<std::size_t>(it - lnE_.begin());
    if (j == 0) j = 1;
    if (j >= lnE_.size()) j = lnE_.size() - 1;

    const std::size_t i = j - 1;

    // Interpolacao linear em (ln E, ln sigma). Uma lei de potencia e
    // uma reta neste plano, entao e reproduzida EXATAMENTE -- e o que
    // o teste T6 verifica.
    const double t = (x - lnE_[i])/(lnE_[j] - lnE_[i]);
    return std::exp(lnSigma_[i] + t*(lnSigma_[j] - lnSigma_[i]));
}

// ---------------------------------------------------------------------

const std::vector<std::string>& required_metadata_keys()
{
    static const std::vector<std::string> k = {
        "convention_y", "target", "projectile", "current",
        "units_sigma", "units_E", "M_Z_GeV", "dipole_model", "generated_by"
    };
    return k;
}

bool TableCrossSection::has_meta(const std::string& chave) const
{
    for (const auto& kv : meta_) if (kv.first == chave) return true;
    return false;
}

const std::string& TableCrossSection::meta(const std::string& chave) const
{
    for (const auto& kv : meta_) if (kv.first == chave) return kv.second;
    throw std::runtime_error("TableCrossSection: metadado ausente '" + chave
                             + "' em " + path_);
}

void TableCrossSection::require_full_metadata() const
{
    for (const std::string& ch : required_metadata_keys()) {
        if (!has_meta(ch)) {
            throw std::runtime_error(
                "TableCrossSection: falta o metadado obrigatorio '" + ch
                + "' em " + path_ + ". Uma tabela sem convencao declarada nao "
                "pode ser comparada com outra: e ai que mora o fator 2.");
        }
    }
}

void assert_composable(const TableCrossSection& cc, const TableCrossSection& nc)
{
    cc.require_full_metadata();
    nc.require_full_metadata();

    static const char* iguais[] = {
        "target", "projectile", "units_sigma", "units_E"
    };
    for (const char* ch : iguais) {
        if (cc.meta(ch) != nc.meta(ch)) {
            std::ostringstream m;
            m << "assert_composable: '" << ch << "' difere entre as tabelas.\n"
              << "  " << cc.path() << " : " << cc.meta(ch) << "\n"
              << "  " << nc.path() << " : " << nc.meta(ch) << "\n"
              << "Somar as duas nao significa nada.";
            throw std::runtime_error(m.str());
        }
    }

    if (cc.meta("current") != "CC" || nc.meta("current") != "NC") {
        std::ostringstream m;
        m << "assert_composable: esperava uma tabela CC e uma NC, e recebi\n"
          << "  " << cc.path() << " : current = " << cc.meta("current") << "\n"
          << "  " << nc.path() << " : current = " << nc.meta("current") << "\n"
          << "Somar duas do mesmo current conta a mesma coisa duas vezes, e o\n"
          << "resultado e uma secao de choque plausivel e errada por um fator 2.";
        throw std::runtime_error(m.str());
    }
}

void write_metadata_header(std::ostream& os, const std::string& tag,
                           const TableCrossSection& t)
{
    os << "# " << tag << ".path = " << t.path() << "\n";
    os << "# " << tag << ".E_min_GeV = " << t.E_min()
       << "   E_max_GeV = " << t.E_max()
       << "   n_pontos = " << t.size() << "\n";
    for (const auto& kv : t.metadata()) {
        os << "# " << tag << "." << kv.first << " = " << kv.second << "\n";
    }
}

void assert_comparable(const TableCrossSection& a, const TableCrossSection& b)
{
    a.require_full_metadata();
    b.require_full_metadata();

    // O que TEM de coincidir para a razao entre as duas significar algo.
    // dipole_model e generated_by ficam de fora: e exatamente ali que as
    // duas devem diferir.
    static const char* iguais[] = {
        "target", "projectile", "current", "units_sigma", "units_E"
    };
    for (const char* ch : iguais) {
        if (a.meta(ch) != b.meta(ch)) {
            std::ostringstream m;
            m << "assert_comparable: '" << ch << "' difere entre as tabelas.\n"
              << "  " << a.path() << " : " << a.meta(ch) << "\n"
              << "  " << b.path() << " : " << b.meta(ch) << "\n"
              << "Comparar as duas nao mede nada.";
            throw std::runtime_error(m.str());
        }
    }
}

// ---------------------------------------------------------------------

double SumCrossSection::E_domain_min() const
{
    double lo = 0.0;
    for (const auto& p : parts_) lo = std::max(lo, p->E_domain_min());
    return lo;
}

double SumCrossSection::E_domain_max() const
{
    double hi = std::numeric_limits<double>::infinity();
    for (const auto& p : parts_) hi = std::min(hi, p->E_domain_max());
    return hi;
}

std::string ScaledCrossSection::name() const
{
    std::ostringstream m;
    m << k_ << " * " << base_->name();
    return m.str();
}

double SumCrossSection::sigma_tot(double E_GeV) const
{
    double s = 0.0;
    for (const auto& p : parts_) s += p->sigma_tot(E_GeV);
    return s;
}

std::string SumCrossSection::name() const
{
    std::ostringstream m;
    m << "Sum[";
    for (std::size_t i = 0; i < parts_.size(); ++i) {
        if (i) m << " + ";
        m << parts_[i]->name();
    }
    m << "]";
    return m.str();
}

} // namespace phasis
