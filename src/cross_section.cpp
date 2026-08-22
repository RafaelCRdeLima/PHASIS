#include "phasis/cross_section.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace phasis {

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

} // namespace phasis
