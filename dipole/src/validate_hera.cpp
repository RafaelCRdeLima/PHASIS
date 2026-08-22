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

#include <cstring>

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
    const IIMParameters& iim,
    bool useIIM,
    int Nr,
    int Nz
)
{
    double FT = 0.0;
    double FL = 0.0;

    for (QuarkFlavor f : flavors) {
        Parameters wf = makeEMParameters(f, masses);

        if (useIIM) {
            FT += FT_IIM(x, Q2, wf, iim, Nr, Nz);
            FL += FL_IIM(x, Q2, wf, iim, Nr, Nz);
        } else {
            FT += FT_GBW(x, Q2, wf, gbw, Nr, Nz);
            FL += FL_GBW(x, Q2, wf, gbw, Nr, Nz);
        }
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
        "data/hera/hera_nc_ep_920.dat";
    std::string out_file = "data/hera_validation.dat";
    std::string model = "GBW";
    double iim_x0 = -1.0;
    int iim_Nb = -1;

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
        else if (a == "--model")   model = argv[++i];
        else if (a == "--iim-x0")  iim_x0 = std::stod(argv[++i]);
        else if (a == "--iim-Nb")  iim_Nb = std::stoi(argv[++i]);
        else {
            std::cerr << "Argumento desconhecido: " << a << "\n";
            return 2;
        }
    }

    std::vector<QuarkFlavor> flavors = {
        QuarkFlavor::u, QuarkFlavor::d, QuarkFlavor::s
    };
    if (flavor_set == "udsc") flavors.push_back(QuarkFlavor::c);

    GBWParameters gbw;
    IIMParameters iim;
    if (iim_x0 > 0.0) iim.x0 = iim_x0;
    if (iim_Nb > 0)   iim.Nb = iim_Nb;

    const bool useIIM = (model == "IIM" || model == "bCGC");

    // massas default = as do ajuste do modelo escolhido
    QuarkMasses masses;
    if (useIIM) {
        masses.u = masses.d = masses.s = iim.m_light;
        masses.c = iim.m_charm;
    } else {
        masses.u = masses.d = masses.s = gbw.m_light;
        masses.c = gbw.m_charm;
    }
    if (mq > 0.0) {
        masses.u = masses.d = masses.s = mq;
    }

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
        if (useIIM) {
            out << "# bCGC        = gamma_s " << iim.gamma_s << ", N0 " << iim.N0
                << ", x0 " << iim.x0 << ", lambda " << iim.lambda
                << ", B " << iim.BCGC << ", Nb " << iim.Nb << "\n";
        } else {
            out << "# GBW         = sigma0 " << gbw.sigma0_mb << " mb, lambda "
                << gbw.lambda << ", x0 " << gbw.x0 << "\n";
        }
        out << "# quadratura  = Nr " << Nr << ", Nz " << Nz << "\n";
        out << "# npontos     = " << sel.size() << "\n";
        out << "# Q2 x y sigma_red_dado err sigma_red_modelo F2_modelo FL_modelo pull\n";
        out << std::scientific;

        double chi2 = 0.0;
        double sum_ratio = 0.0;
        double flmin = 1.0e30, flmax = -1.0e30;

        // CRITERIO
        //
        // chi2 NAO e usado como aprovacao/reprovacao, e a razao e fisica:
        // os dados HERA I+II combinados (2015) tem erro relativo mediano de
        // 2.5%, enquanto GBW e bCGC sao ajustes de 3 a 5 parametros feitos
        // a dados dos anos 1990 com erro de 5-10%. A precisao do dado esta
        // MUITO abaixo da acuracia intrinseca do modelo, entao o chi2 mede
        // as limitacoes conhecidas do modelo (falta de evolucao DGLAP em
        // Q^2 alto), nao se o codigo o implementa certo. Cobrar chi2/dof ~ 1
        // seria exigir que um modelo de 1999 supere a qualidade do proprio
        // ajuste que o definiu.
        //
        // O que testa a IMPLEMENTACAO e se o modelo cai na magnitude e na
        // forma certas:
        //
        //   desvio mediano |mod/dado - 1| <= 15%
        //   <mod/dado> em [0.85, 1.15]
        //   F_L/F_2 em [0.05, 0.25]  (sentinela da quadratura)
        //
        // Isso nao e carimbo: o bCGC com os parametros documentados reprova
        // com folga (razao 0.10). O chi2 continua sendo reportado.
        // Ver CAMPANHA_CORRECAO.md, F1.
        double chi2_core = 0.0;
        int n_core = 0;
        std::vector<double> desvios;

        for (std::size_t i = 0; i < sel.size(); ++i) {
            const HeraPoint& p = sel[i];

            StructureEM F = computeEM(p.x, p.Q2, flavors, masses, gbw, iim, useIIM, Nr, Nz);

            const double yy = p.y * p.y / (1.0 + (1.0 - p.y) * (1.0 - p.y));
            const double model = F.F2 - yy * F.FL;

            const double pull = (model - p.sigma_red) / p.err;
            chi2 += pull * pull;
            sum_ratio += model / p.sigma_red;

            if (p.Q2 >= 1.0 && p.Q2 <= 10.0) {
                chi2_core += pull * pull;
                ++n_core;
            }
            desvios.push_back(std::fabs(model / p.sigma_red - 1.0));

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

        const double chi2_core_ndf = n_core > 0 ? chi2_core/n_core : 1.0e30;

        std::sort(desvios.begin(), desvios.end());
        const double desvio_mediano = desvios[desvios.size()/2];

        const bool passou = (desvio_mediano <= 0.15)
                         && (mean_ratio >= 0.85 && mean_ratio <= 1.15)
                         && (flmin >= 0.05 && flmax <= 0.25);

        std::cout << "\n=== F1: limite EM vs sigma_red do HERA ===\n";
        std::cout << "  arquivo        : " << data_file << "\n";
        std::cout << "  janela         : x < " << xmax
                  << ",  " << Q2min << " <= Q2 <= " << Q2max << " GeV^2\n";
        std::cout << "  pontos         : " << sel.size() << "\n";
        std::cout << "  sabores        : " << flavor_set
                  << "   m_uds = " << masses.u << " GeV\n";
        std::cout << "  modelo         : " << (useIIM ? "bCGC/IIM" : "GBW") << "\n";
        std::cout << "  quadratura     : Nr = " << Nr << ", Nz = " << Nz << "\n";
        std::cout << "  ----------------------------------------\n";
        std::cout << "  chi2/ponto (todos)      : " << chi2_ndf << "\n";
        std::cout << "  chi2/ponto (1<=Q2<=10)  : " << chi2_core_ndf
                  << "   [" << n_core << " pontos]\n";
        std::cout << "  desvio mediano          : " << 100.0*desvio_mediano << " %\n";
        std::cout << "  <modelo/dado>           : " << mean_ratio << "\n";
        std::cout << "  F_L/F2         : " << flmin << " a " << flmax
                  << "   (fisico: 0.05 a 0.25)\n";
        std::cout << "  ----------------------------------------\n";
        std::cout << "  criterio: desvio mediano <= 15%, <mod/dado> em [0.85,1.15],\n";
        std::cout << "            F_L/F2 em [0.05,0.25]. chi2 e reportado, nao cobrado\n";
        std::cout << "            (erro do dado 2.5% << acuracia do modelo ~10%).\n";
        std::cout << "  RESULTADO      : " << (passou ? "OK" : "FALHA") << "\n";
        std::cout << "  saida          : " << out_file << "\n\n";

        return passou ? 0 : 1;
    }
    catch (const std::exception& e) {
        std::cerr << "Erro em validate_hera: " << e.what() << "\n";
        return 2;
    }
}
