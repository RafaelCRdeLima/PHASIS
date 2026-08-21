#include "weak_structure_functions.hpp"

#include <stdexcept>

namespace weak {

WeakStructureFunctions::WeakStructureFunctions(
    const std::string& pdf_set,
    int member
)
{
    pdf_.reset(LHAPDF::mkPDF(pdf_set, member));
}

double WeakStructureFunctions::xf(int pid, double x, double Q2) const
{
    if (x <= 0.0 || x >= 1.0 || Q2 <= 0.0) {
        return 0.0;
    }

    return pdf_->xfxQ2(pid, x, Q2);
}

double WeakStructureFunctions::xF3_CC_isoscalar(
    double x,
    double Q2,
    BeamType beam
) const
{
    const double xu    = xf( 2, x, Q2);
    const double xd    = xf( 1, x, Q2);
    const double xs    = xf( 3, x, Q2);
    const double xc    = xf( 4, x, Q2);

    const double xubar = xf(-2, x, Q2);
    const double xdbar = xf(-1, x, Q2);
    const double xsbar = xf(-3, x, Q2);
    const double xcbar = xf(-4, x, Q2);

    if (beam == BeamType::Neutrino) {
        return xu + xd + 2.0*xs - xubar - xdbar - 2.0*xcbar;
    }

    if (beam == BeamType::AntiNeutrino) {
        return xu + xd + 2.0*xc - xubar - xdbar - 2.0*xsbar;
    }

    return 0.0;
}

BeamType parseBeamType(const std::string& name)
{
    if (name == "nu" || name == "neutrino") {
        return BeamType::Neutrino;
    }

    if (name == "nubar" || name == "antineutrino") {
        return BeamType::AntiNeutrino;
    }

    throw std::runtime_error("Beam inválido. Use --beam nu ou --beam nubar.");
}

}