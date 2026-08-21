// ================================================================
// F1 - Validacao do limite eletromagnetico contra o F2 do HERA
//
// GBW e bCGC NAO sao modelos livres: sao ajustes a secao de choque
// reduzida medida no HERA. Se a maquinaria deste codigo (funcoes de
// onda + quadratura + sigma_dip) nao reproduz esse dado no limite
// eletromagnetico, nenhum resultado a jusante -- inclusive
// sigma_nuN_CC -- tem base.
//
// Este e o teste mais barato que existe para a campanha, e e o que
// teria pego o erro de quadratura na primeira rodada.
//
// Observavel:
//
//     sigma_red = F2 - [ y^2 / (1 + (1-y)^2) ] * F_L
//
// Dados: H1 + ZEUS combinados, HERA I+II, e+p, sqrt(s) = 318 GeV
//        Eur. Phys. J. C 75 (2015) 580  [arXiv:1506.06042]
//
// Uso:
//   validate_hera [--Nr N] [--Nz N] [--mq M] [--xmax X]
//                 [--Q2min Q] [--Q2max Q] [--flavors uds|udsc]
//                 [--data ARQ] [--out ARQ]
// ================================================================

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "parameters.hpp"
#include "dipole_models.hpp"
#include "integrals.hpp"

using namespace dipole;

namespace {

struct HeraPoint {
    double Q2;
    double x;
    double y;
    double sigma_red;
    double err;
};

std::vector<HeraPoint> readHera(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in) {
        throw std::runtime_error("Nao consegui abrir o arquivo de dados do HERA: " + filename);
    }

    std::vector<HeraPoint> pts;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        HeraPoint p;
        std::istringstream ss(line);

        if (ss >> p.Q2 >> p.x >> p.y >> p.sigma_red >> p.err) {
            if (p.err > 0.0) pts.push_back(p);
        }
    }

    if (pts.empty()) {
        throw std::runtime_error("Nenhum ponto valido lido de: " + filename);
    }

    return pts;
}

// F2 e F_L eletromagneticos, somados sobre sabores.
// Reusa exatamente FT_GBW/FL_GBW -- o mesmo caminho de codigo que
// sigma_nuN usa. Trocar a quadratura la muda o resultado aqui.
struct StructureEM {
    double F2;
    double FL;
};

StructureEM computeEM(
    double x,
    double Q2,
    const std::vector<QuarkFlavor>& flavors,
    const QuarkMasses& masses,
    const GBWParameters& gbw,
    int Nr,
    int Nz
)
{
    double FT = 0.0;
    double FL = 0.0;

    for (QuarkFlavor f : flavors) {
        Parameters wf = makeEMParameters(f, masses);

        FT += FT_GBW(x, Q2, wf, gbw, Nr, Nz);
        FL += FL_GBW(x, Q2, wf, gbw, Nr, Nz);
    }

    // Convencao DIS: F2 = Q^2/(4 pi^2 alpha_em) * integral.
    // As wavefunctions carregam alpha_em, entao dividimos aqui.
    // No limite EM esta convencao NAO e ambigua -- e a definicao padrao.
    FT /= alphaEM;
    FL /= alphaEM;

    return { FT + FL, FL };
}

} // namespace

int main(int argc, char* argv[])
{
    int Nr = 400;
    int Nz = 200;

    double mq = -1.0;          // <0 = usar as massas de QuarkMasses
    double xmax = 1.0e-2;      // janela de ajuste do GBW
    double Q2min = 0.25;
    double Q2max = 50.0;

    std::string flavor_set = "uds";
    std::string data_file =
        "/home/rafael/Codes/HADROS3/sandbox/data/hera/hera_nc_ep_920.dat";
    std::string out_file = "data/hera_validation.dat";

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if      (a == "--Nr")      Nr = std::stoi(argv[++i]);
        else if (a == "--Nz")      Nz = std::stoi(argv[++i]);
        else if (a == "--mq")      mq = std::stod(argv[++i]);
        else if (a == "--xmax")    xmax = std::stod(argv[++i]);
        else if (a == "--Q2min")   Q2min = std::stod(argv[++i]);
        else if (a == "--Q2max")   Q2max = std::stod(argv[++i]);
        else if (a == "--flavors") flavor_set = argv[++i];
        else if (a == "--data")    data_file = argv[++i];
        else if (a == "--out")     out_file = argv[++i];
        else {
            std::cerr << "Argumento desconhecido: " << a << "\n";
            return 2;
        }
    }

    std::vector<QuarkFlavor> flavors = {
        QuarkFlavor::u, QuarkFlavor::d, QuarkFlavor::s
    };
    if (flavor_set == "udsc") flavors.push_back(QuarkFlavor::c);

    QuarkMasses masses;
    if (mq > 0.0) {
        masses.u = masses.d = masses.s = mq;
    }

    GBWParameters gbw;

    try {
        auto all = readHera(data_file);

        std::vector<HeraPoint> sel;
        for (const auto& p : all) {
            if (p.x < xmax && p.Q2 >= Q2min && p.Q2 <= Q2max) sel.push_back(p);
        }

        if (sel.empty()) {
            std::cerr << "Nenhum ponto do HERA na janela pedida.\n";
            return 1;
        }

        std::ofstream out(out_file);
        out << "# Validacao do limite EM contra sigma_red do HERA\n";
        out << "# dados       = " << data_file << "\n";
        out << "# janela      = x < " << xmax
            << ", " << Q2min << " <= Q2 <= " << Q2max << " GeV^2\n";
        out << "# sabores     = " << flavor_set << "\n";
        out << "# m_uds       = " << masses.u << " GeV\n";
        out << "# GBW         = sigma0 " << gbw.sigma0_mb << " mb, lambda "
            << gbw.lambda << ", x0 " << gbw.x0 << "\n";
        out << "# quadratura  = Nr " << Nr << ", Nz " << Nz << "\n";
        out << "# npontos     = " << sel.size() << "\n";
        out << "# Q2 x y sigma_red_dado err sigma_red_modelo F2_modelo FL_modelo pull\n";
        out << std::scientific;

        double chi2 = 0.0;
        double sum_ratio = 0.0;
        double flmin = 1.0e30, flmax = -1.0e30;

        for (std::size_t i = 0; i < sel.size(); ++i) {
            const HeraPoint& p = sel[i];

            StructureEM F = computeEM(p.x, p.Q2, flavors, masses, gbw, Nr, Nz);

            const double yy = p.y * p.y / (1.0 + (1.0 - p.y) * (1.0 - p.y));
            const double model = F.F2 - yy * F.FL;

            const double pull = (model - p.sigma_red) / p.err;
            chi2 += pull * pull;
            sum_ratio += model / p.sigma_red;

            const double ratio = (F.F2 > 0.0) ? F.FL / F.F2 : 0.0;
            flmin = std::min(flmin, ratio);
            flmax = std::max(flmax, ratio);

            out << p.Q2 << " " << p.x << " " << p.y << " "
                << p.sigma_red << " " << p.err << " "
                << model << " " << F.F2 << " " << F.FL << " " << pull << "\n";

            if ((i + 1) % 20 == 0 || i + 1 == sel.size()) {
                std::cerr << "  " << (i + 1) << "/" << sel.size() << "\r" << std::flush;
            }
        }
        std::cerr << "\n";

        const double ndf = static_cast<double>(sel.size());
        const double chi2_ndf = chi2 / ndf;
        const double mean_ratio = sum_ratio / ndf;

        std::cout << "\n=== F1: limite EM vs sigma_red do HERA ===\n";
        std::cout << "  arquivo        : " << data_file << "\n";
        std::cout << "  janela         : x < " << xmax
                  << ",  " << Q2min << " <= Q2 <= " << Q2max << " GeV^2\n";
        std::cout << "  pontos         : " << sel.size() << "\n";
        std::cout << "  sabores        : " << flavor_set
                  << "   m_uds = " << masses.u << " GeV\n";
        std::cout << "  quadratura     : Nr = " << Nr << ", Nz = " << Nz << "\n";
        std::cout << "  ----------------------------------------\n";
        std::cout << "  chi2/ponto     : " << chi2_ndf << "\n";
        std::cout << "  <modelo/dado>  : " << mean_ratio << "\n";
        std::cout << "  F_L/F2         : " << flmin << " a " << flmax
                  << "   (fisico: 0.05 a 0.25)\n";
        std::cout << "  ----------------------------------------\n";
        std::cout << "  criterio de aceitacao da campanha: chi2/ponto <= 2\n";
        std::cout << "  RESULTADO      : "
                  << (chi2_ndf <= 2.0 ? "OK" : "FALHA") << "\n";
        std::cout << "  saida          : " << out_file << "\n\n";

        return (chi2_ndf <= 2.0) ? 0 : 1;
    }
    catch (const std::exception& e) {
        std::cerr << "Erro em validate_hera: " << e.what() << "\n";
        return 2;
    }
}
