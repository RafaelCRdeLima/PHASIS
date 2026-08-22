#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <functional>
#include <stdexcept>

#include "parameters.hpp"
#include "dipole_models.hpp"
#include "integrals.hpp"
#include "weak_structure_functions.hpp"
#include <memory>

using namespace dipole;

// Conversão:
// 1 GeV^{-2} = 0.389379 mb = 0.389379e-27 cm^2
constexpr double GeVminus2_to_cm2 = 0.389379e-27;
constexpr double MN = 0.938272; // massa do núcleon em GeV

enum class DipoleModel {
    GBW,
    IIM
};

struct StructureFunctions {
    double FT = 0.0;
    double FL = 0.0;
    double F2 = 0.0;
};

double largeXFactor(double x)
{
    if (x >= 1.0) return 0.0;
    if (x <= 0.0) return 0.0;
    return std::pow(1.0 - x, 7.0);
}

DipoleModel parseDipoleModel(const std::string& model_name)
{
    if (model_name == "GBW") return DipoleModel::GBW;
    if (model_name == "IIM") return DipoleModel::IIM;

    throw std::runtime_error("Modelo inválido. Use --model GBW ou --model IIM.");
}

StructureFunctions computeChannelCC(
    double x,
    double Q2,
    QuarkFlavor q,
    QuarkFlavor qbar,
    const GBWParameters& gbw,
    const IIMParameters& iim,
    DipoleModel model,
    int Nr,
    int Nz
)
{
    Parameters wf = makeCCParameters(q, qbar);

    double FT = 0.0;
    double FL = 0.0;

    if (model == DipoleModel::GBW) {
        FT = FT_GBW(x, Q2, wf, gbw, Nr, Nz);
        FL = FL_GBW(x, Q2, wf, gbw, Nr, Nz);
    }
    else if (model == DipoleModel::IIM) {
        FT = FT_IIM(x, Q2, wf, iim, Nr, Nz);
        FL = FL_IIM(x, Q2, wf, iim, Nr, Nz);
    }

    // Convenção DIS: remove alphaEW global das wavefunctions.
    FT /= wf.alphaEW;
    FL /= wf.alphaEW;

    StructureFunctions F;
    F.FT = FT;
    F.FL = FL;
    F.F2 = FT + FL;

    return F;
}

StructureFunctions computeFCCtotal(
    double x,
    double Q2,
    const GBWParameters& gbw,
    const IIMParameters& iim,
    DipoleModel model,
    int Nr,
    int Nz
)
{
    StructureFunctions ud = computeChannelCC(
        x, Q2,
        QuarkFlavor::u,
        QuarkFlavor::d,
        gbw,
        iim,
        model,
        Nr,
        Nz
    );

    StructureFunctions cs = computeChannelCC(
        x, Q2,
        QuarkFlavor::c,
        QuarkFlavor::s,
        gbw,
        iim,
        model,
        Nr,
        Nz
    );

    const double lx = largeXFactor(x);

    StructureFunctions total;
    total.FT = lx*(ud.FT + cs.FT);
    total.FL = lx*(ud.FL + cs.FL);
    total.F2 = lx*(ud.F2 + cs.F2);

    return total;
}

double d2sigma_dxdy_CC(
    double Enu,
    double x,
    double Q2,
    const GBWParameters& gbw,
    const IIMParameters& iim,
    DipoleModel model,
    int Nr,
    int Nz,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam
)
{
    const double s = 2.0*MN*Enu;
    const double y = Q2/(x*s);

    if (x <= 0.0 || x >= 1.0) return 0.0;
    if (y <= 0.0 || y >= 1.0) return 0.0;

    StructureFunctions F =
        computeFCCtotal(x, Q2, gbw, iim, model, Nr, Nz);

    const double propagator =
        std::pow(MW*MW/(Q2 + MW*MW), 2);

    const double prefactor =
        GF*GF*MN*Enu/pi * propagator;

    double xF3 = 0.0;

    if (useF3 && weakSF != nullptr) {
        xF3 = weakSF->xF3_CC_isoscalar(x, Q2, beam);
    }

    const double f3_sign =
        (beam == weak::BeamType::Neutrino) ? 1.0 : -1.0;

    const double bracket =
        0.5*(1.0 + std::pow(1.0-y, 2))*F.F2
        - 0.5*y*y*F.FL
        + f3_sign*y*(1.0 - 0.5*y)*xF3;

    return prefactor*bracket;
}

double sigmaNuN_CC(
    double Enu,
    const GBWParameters& gbw,
    const IIMParameters& iim,
    DipoleModel model,
    int NlogQ,
    int Nlogx,
    int Nr,
    int Nz,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam
)
{
    const double s = 2.0*MN*Enu;

    const double Q2min = 1.0;
    const double Q2max = 0.999*s;

    if (Q2max <= Q2min) return 0.0;

    const double logQmin = std::log(Q2min);
    const double logQmax = std::log(Q2max);

    auto integrandLogQ = [&](double logQ2)
    {
        const double Q2 = std::exp(logQ2);

        const double xmin = Q2/s;

        double xmax = 0.999999;

        if (model == DipoleModel::IIM) {
            xmax = 1.0e-2;
        }

        if (xmin >= xmax) return 0.0;

        const double logxmin = std::log(xmin);
        const double logxmax = std::log(xmax);

        auto integrandLogx = [&](double logx)
        {
            const double x = std::exp(logx);

            const double dsigma =
                     d2sigma_dxdy_CC(
                         Enu,
                         x,
                         Q2,
                         gbw,
                         iim,
                         model,
                         Nr,
                         Nz,
                         useF3,
                         weakSF,
                         beam
                     );

            return dsigma/s;
        };

        const double inner =
            simpson(
                integrandLogx,
                logxmin,
                logxmax,
                Nlogx
            );

        return Q2*inner;
    };

    return simpson(
        integrandLogQ,
        logQmin,
        logQmax,
        NlogQ
    );
}

int main(int argc, char* argv[])
{
    bool useF3 = false;
    std::string pdf_set = "CT10nlo"; //"NNPDF31_nlo_as_0118";
    std::string beam_name = "nu";
    double logEmin = 3.0;
    double logEmax = 14.0;
    int NE = 20;

    int NlogQ = 16;
    int Nlogx = 16;

    int Nr = 200;
    int Nz = 200;

    std::string model_name = "GBW";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--model") model_name = argv[++i];

        else if (arg == "--logEmin") logEmin = std::stod(argv[++i]);
        else if (arg == "--logEmax") logEmax = std::stod(argv[++i]);
        else if (arg == "--NE") NE = std::stoi(argv[++i]);

        else if (arg == "--NlogQ") NlogQ = std::stoi(argv[++i]);
        else if (arg == "--Nlogx") Nlogx = std::stoi(argv[++i]);

        else if (arg == "--Nr") Nr = std::stoi(argv[++i]);
        else if (arg == "--Nz") Nz = std::stoi(argv[++i]);
        else if (arg == "--use-F3") useF3 = std::stoi(argv[++i]) != 0;
        else if (arg == "--pdf-set") pdf_set = argv[++i];
        else if (arg == "--beam") beam_name = argv[++i];
    }

    DipoleModel model = parseDipoleModel(model_name);

    GBWParameters gbw;
    IIMParameters iim;

    weak::BeamType beam = weak::parseBeamType(beam_name);

    std::unique_ptr<weak::WeakStructureFunctions> weakSF;

    if (useF3) {
        weakSF = std::make_unique<weak::WeakStructureFunctions>(pdf_set, 0);
    }

    std::string output_file =
        "data/sigma_nuN_CC_" + model_name + ".dat";

    std::ofstream out(output_file);

    out << "# model " << model_name << "\n";
    out << "# Enu_GeV sigma_GeV_minus2 sigma_cm2\n";

    for (int i = 0; i < NE; ++i) {

        const double logE =
            logEmin + i*(logEmax - logEmin)/(NE - 1);

        const double Enu = std::pow(10.0, logE);

        const double sigma_gev2 =
        sigmaNuN_CC(
            Enu,
            gbw,
            iim,
            model,
            NlogQ,
            Nlogx,
            Nr,
            Nz,
            useF3,
            weakSF.get(),
            beam
        );

        const double sigma_cm2 =
            sigma_gev2 * GeVminus2_to_cm2;

        out << Enu << " "
            << sigma_gev2 << " "
            << sigma_cm2 << "\n";

        std::cout << "model = " << model_name
                  << "   E_nu = " << Enu
                  << " GeV   sigma_CC = "
                  << sigma_cm2
                  << " cm^2\n";
    }

    std::cout << "Arquivo gerado: " << output_file << "\n";

    return 0;
}