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

double WeakStructureFunctions::xF3_NC_isoscalar(
    double x,
    double Q2,
    BeamType beam
) const
{
    // g_V g_A por tipo, com sin^2(thetaW) = 0.23122 (parameters.hpp).
    // Escritos aqui como literais para nao arrastar parameters.hpp
    // para dentro do wrapper do LHAPDF.
    constexpr double s2w = 0.23122;
    const double gVu =  0.5 - (4.0/3.0)*s2w, gAu =  0.5;
    const double gVd = -0.5 + (2.0/3.0)*s2w, gAd = -0.5;
    const double Gu = gVu*gAu;
    const double Gd = gVd*gAd;

    const double xu = impl_->xf( 2, x, Q2), xub = impl_->xf(-2, x, Q2);
    const double xd = impl_->xf( 1, x, Q2), xdb = impl_->xf(-1, x, Q2);
    const double xs = impl_->xf( 3, x, Q2), xsb = impl_->xf(-3, x, Q2);
    const double xc = impl_->xf( 4, x, Q2), xcb = impl_->xf(-4, x, Q2);

    // Alvo isoescalar: u_N = d_N = (u_p + d_p)/2, logo a valencia de
    // u_N e a de d_N sao ambas V/2.
    const double V = (xu + xd) - (xub + xdb);

    const double xF3 = (Gu + Gd)*V + 2.0*Gd*(xs - xsb) + 2.0*Gu*(xc - xcb);

    return (beam == BeamType::Neutrino) ? xF3 : -xF3;
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

double WeakStructureFunctions::xF3_NC_isoscalar(double, double, BeamType) const
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
