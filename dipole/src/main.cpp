#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "parameters.hpp"
#include "wavefunctions.hpp"
#include "dipole_models.hpp"
#include "integrals.hpp"

using namespace dipole;

QuarkFlavor parseFlavor(const std::string& s)
{
    if (s == "u") return QuarkFlavor::u;
    if (s == "d") return QuarkFlavor::d;
    if (s == "s") return QuarkFlavor::s;
    if (s == "c") return QuarkFlavor::c;
    if (s == "b") return QuarkFlavor::b;
    if (s == "t") return QuarkFlavor::t;
    throw std::runtime_error("Flavor inválido.");
}

int main(int argc, char* argv[])
{
    std::string mode = "wavefunctions";
    std::string current = "NC";
    std::string flavor = "u";

    double Q2 = 10.0;
    double x  = 1.0e-5;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--mode") mode = argv[++i];
        else if (arg == "--current") current = argv[++i];
        else if (arg == "--flavor") flavor = argv[++i];
        else if (arg == "--Q2") Q2 = std::stod(argv[++i]);
        else if (arg == "--x") x = std::stod(argv[++i]);
    }

    QuarkFlavor q = parseFlavor(flavor);

    Parameters wf;

    if (current == "NC") {
        wf = makeNCParameters(q);
    } else if (current == "CC") {
        // exemplo padrão: W -> c sbar
        wf = makeCCParameters(QuarkFlavor::c, QuarkFlavor::s);
    } else {
        throw std::runtime_error("Current deve ser NC ou CC.");
    }

    GBWParameters gbw;

    if (mode == "wavefunctions") {

        std::ofstream out("data/wavefunctions.dat");

        out << "# r z psiT2 psiL2\n";

        for (int iz = 0; iz <= 100; ++iz) {
            double z = wf.zMin + iz*(wf.zMax - wf.zMin)/100.0;

            for (int ir = 0; ir <= 300; ++ir) {
                double log_r = std::log10(wf.rMin)
                             + ir*(std::log10(wf.rMax) - std::log10(wf.rMin))/300.0;

                double r = std::pow(10.0, log_r);

                out << r << " "
                    << z << " "
                    << psiT2(r,z,Q2,wf) << " "
                    << psiL2(r,z,Q2,wf) << "\n";
            }

            out << "\n";
        }

        std::cout << "Arquivo gerado: data/wavefunctions.dat\n";
    }

    else if (mode == "integrand") {

        std::ofstream out("data/integrand.dat");

        out << "# r z integrandT integrandL\n";

        for (int iz = 0; iz <= 100; ++iz) {
            double z = wf.zMin + iz*(wf.zMax - wf.zMin)/100.0;

            for (int ir = 0; ir <= 300; ++ir) {
                double log_r = std::log10(wf.rMin)
                             + ir*(std::log10(wf.rMax) - std::log10(wf.rMin))/300.0;

                double r = std::pow(10.0, log_r);

                out << r << " "
                    << z << " "
                    << integrandT(r,z,x,Q2,wf,gbw) << " "
                    << integrandL(r,z,x,Q2,wf,gbw) << "\n";
            }

            out << "\n";
        }

        std::cout << "Arquivo gerado: data/integrand.dat\n";
    }

    else if (mode == "sigma") {

        double FT = FT_GBW(x,Q2,wf,gbw);
        double FL = FL_GBW(x,Q2,wf,gbw);
        double F2 = FT + FL;

        std::ofstream out("data/structure_functions.dat");

        out << "# x Q2 FT FL F2\n";
        out << x << " " << Q2 << " "
            << FT << " " << FL << " " << F2 << "\n";

        std::cout << "FT = " << FT << "\n";
        std::cout << "FL = " << FL << "\n";
        std::cout << "F2 = " << F2 << "\n";
        std::cout << "Arquivo gerado: data/structure_functions.dat\n";
    }

    else {
        throw std::runtime_error("Modo inválido.");
    }

    return 0;
}