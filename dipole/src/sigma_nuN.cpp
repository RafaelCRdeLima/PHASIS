#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "parameters.hpp"
#include "dipole_models.hpp"
#include "integrals.hpp"
#include "structure_table.hpp"
#include "sigma_nuN_core.hpp"
#include "weak_structure_functions.hpp"

using namespace dipole;

enum class DipoleModel { GBW, IIM };

static DipoleModel parseDipoleModel(const std::string& model_name)
{
    if (model_name == "GBW") return DipoleModel::GBW;
    if (model_name == "IIM") return DipoleModel::IIM;
    throw std::runtime_error("Modelo invalido. Use --model GBW ou --model IIM.");
}

int main(int argc, char* argv[])
{
    bool useF3 = false;
    std::string pdf_set = "CT10nlo";
    std::string beam_name = "nu";
    double logEmin = 3.0;
    double logEmax = 14.0;
    int NE = 20;

    // F4: densidade de nos por decada, em vez de numero fixo.
    double nodesQ = 16.0;
    double nodesX = 16.0;

    // Quadratura interna em (r,z) usada para montar a tabela.
    QuadratureGrid quad;
    quad.Nr = 120;
    quad.Nz = 120;

    TableSpec spec;

    double Q2min = 1.0;
    std::string model_name = "GBW";
    bool verbose = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if      (arg == "--model")    model_name = argv[++i];
        else if (arg == "--logEmin")  logEmin = std::stod(argv[++i]);
        else if (arg == "--logEmax")  logEmax = std::stod(argv[++i]);
        else if (arg == "--NE")       NE = std::stoi(argv[++i]);
        else if (arg == "--nodesQ")   nodesQ = std::stod(argv[++i]);
        else if (arg == "--nodesX")   nodesX = std::stod(argv[++i]);
        else if (arg == "--Nr")       quad.Nr = std::stoi(argv[++i]);
        else if (arg == "--Nz")       quad.Nz = std::stoi(argv[++i]);
        else if (arg == "--Ntabx")    spec.Nx = std::stoi(argv[++i]);
        else if (arg == "--NtabQ")    spec.NQ = std::stoi(argv[++i]);
        else if (arg == "--Q2min")    Q2min = std::stod(argv[++i]);
        else if (arg == "--use-F3")   useF3 = std::stoi(argv[++i]) != 0;
        else if (arg == "--pdf-set")  pdf_set = argv[++i];
        else if (arg == "--beam")     beam_name = argv[++i];
        else if (arg == "--quiet")    verbose = false;
        else {
            std::cerr << "Argumento desconhecido: " << arg << "\n";
            return 2;
        }
    }

    if (NE < 1) {
        std::cerr << "--NE precisa ser >= 1.\n";
        return 2;
    }

    DipoleModel model = parseDipoleModel(model_name);

    GBWParameters gbw;
    IIMParameters iim;
    QuarkMasses masses;

    weak::BeamType beam = weak::parseBeamType(beam_name);

    std::unique_ptr<weak::WeakStructureFunctions> weakSF;
    if (useF3) {
        weakSF = std::make_unique<weak::WeakStructureFunctions>(pdf_set, 0);
    }

    // A tabela precisa cobrir todo o (x, Q^2) que o laco de energia
    // vai pedir: x_min = Q2min/s_max e Q2_max = 0.999*s_max.
    const double s_max = 2.0*MN*std::pow(10.0, logEmax);
    spec.Q2Min = Q2min;
    spec.Q2Max = 1.05*0.999*s_max;
    spec.xMin  = 0.5*Q2min/s_max;
    spec.xMax  = 1.0;

    // Kutak-Kwiecinski: os dipolos favorecidos por Cabibbo sao
    // ud~(du~) e cs~(sc~). Os dois canais somados dao o Sigma(carga^2)
    // efetivo = 4 das eqs. (12) e (13).
    const std::vector<StructureTable::Channel> channels = {
        { QuarkFlavor::u, QuarkFlavor::d },
        { QuarkFlavor::c, QuarkFlavor::s }
    };

    if (verbose) {
        std::cerr << "Montando a tabela de F_T, F_L (uma vez, serve todas as energias)\n";
    }

    const StructureTable table(
        model == DipoleModel::GBW ? DipoleModelId::GBW : DipoleModelId::IIM,
        channels, gbw, iim, quad, spec, masses, verbose
    );

    if (verbose) {
        std::cerr << "  " << table.describe() << "\n"
                  << "  construida em " << table.buildSeconds() << " s\n\n";
    }

    const std::string output_file =
        "data/sigma_nuN_CC_" + model_name + ".dat";

    std::ofstream out(output_file);
    out << std::setprecision(10);

    // F9: proveniencia. A tabela tem de poder reproduzir a si mesma
    // so pelo cabecalho.
    out << "# model " << model_name << "\n";
    out << "# gerado_por sigma_nuN (PHASIS/dipole)\n";
    out << "# formalismo Kutak-Kwiecinski EPJ C29 (2003) 521, eqs (2),(8),(9),(12),(13)\n";
    out << "# beam " << beam_name << "\n";
    out << "# use_F3 " << (useF3 ? 1 : 0);
    if (useF3) out << "  pdf_set " << pdf_set;
    out << "\n";
    out << "# Q2min_GeV2 " << Q2min << "\n";
    out << "# largeXFactor (1-x)^7   [regra de contagem de constituintes, n_s=4]\n";
    if (model == DipoleModel::GBW) {
        out << "# GBW sigma0_mb " << gbw.sigma0_mb
            << "  lambda " << gbw.lambda
            << "  x0 " << gbw.x0
            << "  Q0sq " << gbw.Q0sq << "\n";
    } else {
        out << "# bCGC B_CGC " << iim.BCGC
            << "  gamma_s " << iim.gamma_s
            << "  N0 " << iim.N0
            << "  x0 " << iim.x0
            << "  lambda " << iim.lambda
            << "  bMax " << iim.bMax
            << "  Nb " << iim.Nb << "\n";
    }
    out << "# massas_GeV  u " << masses.u << "  d " << masses.d
        << "  s " << masses.s << "  c " << masses.c << "\n";
    out << "# canais ud, cs\n";
    out << "# quadratura_rz  Nr " << quad.Nr << "  Nz " << quad.Nz
        << "  rMin " << quad.rMin << "  rMax " << quad.rMax
        << "  zMin " << quad.zMin << "   [ln r; z log nas duas pontas]\n";
    out << "# tabela_F  Nx " << spec.Nx << "  NQ " << spec.NQ
        << "  xMin " << spec.xMin << "  xMax " << spec.xMax
        << "  Q2Min " << spec.Q2Min << "  Q2Max " << spec.Q2Max << "\n";
    out << "# integracao_sigma  nos_por_decada_Q " << nodesQ
        << "  nos_por_decada_x " << nodesX << "\n";
    out << "# Enu_GeV sigma_GeV_minus2 sigma_cm2\n";

    for (int i = 0; i < NE; ++i) {

        const double logE = (NE == 1)
            ? logEmin
            : logEmin + i*(logEmax - logEmin)/(NE - 1);

        const double Enu = std::pow(10.0, logE);

        int nodesQused = 0;
        const double sigma_gev2 = sigmaNuN_CC(
            Enu, table, nodesQ, nodesX, useF3, weakSF.get(), beam,
            Q2min, &nodesQused
        );

        const double sigma_cm2 = sigma_gev2*GeVminus2_to_cm2;

        out << Enu << " " << sigma_gev2 << " " << sigma_cm2 << "\n";

        if (verbose) {
            std::cout << "model = " << model_name
                      << "   E_nu = " << Enu
                      << " GeV   sigma_CC = " << sigma_cm2
                      << " cm^2   (NlogQ = " << nodesQused << ")\n";
        }
    }

    if (table.clampedQueries() > 0) {
        std::cerr << "AVISO: " << table.clampedQueries()
                  << " consultas fora da faixa da tabela foram truncadas.\n";
    }

    std::cout << "Arquivo gerado: " << output_file << "\n";
    return 0;
}
