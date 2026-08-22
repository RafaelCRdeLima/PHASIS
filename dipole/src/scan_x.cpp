#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <memory>

#include "parameters.hpp"
#include "dipole_models.hpp"
#include "integrals.hpp"
#include "weak_structure_functions.hpp"

using namespace dipole;

int main(int argc, char* argv[])
{
    bool useF3 = true;
    std::string pdf_set = "NNPDF31_nlo_as_0118";
    std::string beam_name = "nu";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--use-F3") {
            useF3 = std::stoi(argv[++i]) != 0;
        }
        else if (arg == "--pdf-set") {
            pdf_set = argv[++i];
        }
        else if (arg == "--beam") {
            beam_name = argv[++i];
        }
    }

    Parameters wf = makeCCParameters(QuarkFlavor::c, QuarkFlavor::s);

    GBWParameters gbw;
    IIMParameters iim;

    const double Q2 = MW*MW;

    const double logxmin = -8.0;
    const double logxmax = -2.0;
    const int Nx = 80;

    weak::BeamType beam = weak::parseBeamType(beam_name);

    std::unique_ptr<weak::WeakStructureFunctions> weakSF;

    if (useF3) {
        weakSF = std::make_unique<weak::WeakStructureFunctions>(pdf_set, 0);
    }

    std::ofstream out("data/scan_x_F2_CC_models.dat");

    out << "# x Q2 FT_GBW FL_GBW F2_GBW xF3 FT_IIM FL_IIM F2_IIM\n";

    for (int i = 0; i <= Nx; ++i) {
        double logx = logxmin + i*(logxmax - logxmin)/Nx;
        double x = std::pow(10.0, logx);

        // F5: convencao DIS -- as wavefunctions carregam alphaEW e as
        // funcoes de estrutura de Kutak-Kwiecinski (eq. 9) nao. Sem esta
        // divisao, o F2 escrito aqui saia 1/alpha_CC = 236x menor que o
        // que entra na secao de choque, e era publicado lado a lado com
        // xF3 do LHAPDF, que esta na normalizacao padrao.
        // Ver CAMPANHA_CORRECAO.md secao 4, Q1.
        double FT_gbw = FT_GBW(x, Q2, wf, gbw, 200, 200)/wf.alphaEW;
        double FL_gbw = FL_GBW(x, Q2, wf, gbw, 200, 200)/wf.alphaEW;
        double F2_gbw = FT_gbw + FL_gbw;

        double FT_iim = FT_IIM(x, Q2, wf, iim, 200, 200)/wf.alphaEW;
        double FL_iim = FL_IIM(x, Q2, wf, iim, 200, 200)/wf.alphaEW;
        double F2_iim = FT_iim + FL_iim;

        double xF3 = 0.0;

        if (useF3 && weakSF) {
            xF3 = weakSF->xF3_CC_isoscalar(x, Q2, beam);
        }

        out << x << " " << Q2 << " "
            << FT_gbw << " " << FL_gbw << " " << F2_gbw << " "
            << xF3 << " "
            << FT_iim << " " << FL_iim << " " << F2_iim << "\n";

        std::cout << "x = " << x
                  << "  F2_GBW = " << F2_gbw
                  << "  xF3 = " << xF3
                  << "  F2_IIM = " << F2_iim
                  << "\n";
    }

    std::cout << "Arquivo gerado: data/scan_x_F2_CC_models.dat\n";

    return 0;
}