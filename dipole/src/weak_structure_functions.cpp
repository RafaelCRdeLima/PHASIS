#include "weak_structure_functions.hpp"

#include <stdexcept>

#ifdef WITH_LHAPDF
#include <LHAPDF/LHAPDF.h>
#endif

namespace weak {

#ifdef WITH_LHAPDF

struct WeakStructureFunctions::Impl {
    std::unique_ptr<LHAPDF::PDF> pdf;

    double xf(int pid, double x, double Q2) const
    {
        if (x <= 0.0 || x >= 1.0 || Q2 <= 0.0) return 0.0;
        return pdf->xfxQ2(pid, x, Q2);
    }
};

bool WeakStructureFunctions::available() { return true; }

WeakStructureFunctions::WeakStructureFunctions(
    const std::string& pdf_set,
    int member
)
    : impl_(new Impl)
{
    impl_->pdf.reset(LHAPDF::mkPDF(pdf_set, member));
}

WeakStructureFunctions::~WeakStructureFunctions() = default;

double WeakStructureFunctions::xF3_CC_isoscalar(
    double x,
    double Q2,
    BeamType beam
) const
{
    const double xu    = impl_->xf( 2, x, Q2);
    const double xd    = impl_->xf( 1, x, Q2);
    const double xs    = impl_->xf( 3, x, Q2);
    const double xc    = impl_->xf( 4, x, Q2);

    const double xubar = impl_->xf(-2, x, Q2);
    const double xdbar = impl_->xf(-1, x, Q2);
    const double xsbar = impl_->xf(-3, x, Q2);
    const double xcbar = impl_->xf(-4, x, Q2);

    if (beam == BeamType::Neutrino) {
        return xu + xd + 2.0*xs - xubar - xdbar - 2.0*xcbar;
    }

    if (beam == BeamType::AntiNeutrino) {
        return xu + xd + 2.0*xc - xubar - xdbar - 2.0*xsbar;
    }

    return 0.0;
}

#else   // sem LHAPDF

struct WeakStructureFunctions::Impl {};

bool WeakStructureFunctions::available() { return false; }

WeakStructureFunctions::WeakStructureFunctions(const std::string&, int)
{
    throw std::runtime_error(
        "Este binario foi compilado SEM LHAPDF, entao --use-F3 1 nao esta "
        "disponivel. Rode com --use-F3 0, ou recompile com o LHAPDF no PATH "
        "(o Makefile detecta lhapdf-config sozinho)."
    );
}

WeakStructureFunctions::~WeakStructureFunctions() = default;

double WeakStructureFunctions::xF3_CC_isoscalar(double, double, BeamType) const
{
    return 0.0;
}

#endif

BeamType parseBeamType(const std::string& name)
{
    if (name == "nu" || name == "neutrino") {
        return BeamType::Neutrino;
    }
    if (name == "nubar" || name == "antineutrino") {
        return BeamType::AntiNeutrino;
    }
    throw std::runtime_error("Beam invalido. Use --beam nu ou --beam nubar.");
}

}
