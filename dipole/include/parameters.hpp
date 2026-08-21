#ifndef DIPOLE_PARAMETERS_HPP
#define DIPOLE_PARAMETERS_HPP

#include <cmath>
#include <stdexcept>
#include <string>

namespace dipole {

constexpr double pi = 3.14159265358979323846;

// Constantes em GeV
constexpr double GF = 1.1663787e-5;   // GeV^{-2}
constexpr double MW = 80.379;         // GeV
constexpr double MZ = 91.1876;        // GeV
constexpr double sin2ThetaW = 0.23122;

// Massas efetivas dos quarks em GeV
struct QuarkMasses {
    double u = 0.03;
    double d = 0.03;
    double s = 0.03;
    double c = 1.3113;
    double b = 4.75;
    double t = 173.0;
};

enum class CurrentType {
    CC,
    NC
};

enum class QuarkFlavor {
    u, d, s, c, b, t
};

struct ElectroweakCouplings {
    double gV;
    double gA;
    double alphaEW;
};

// Parâmetros gerais do cálculo
struct Parameters {
    CurrentType current = CurrentType::NC;
    QuarkFlavor flavor = QuarkFlavor::u;

    // massas dos quarks no dipolo
    double m = 0.03;
    double mu = 0.03;

    // acoplamentos fracos
    double gV = 0.5;
    double gA = 0.5;

    // alpha_ew = g_W^2/(4 pi) ou g_Z^2/(4 pi)
    double alphaEW = (MZ*MZ*GF/std::sqrt(2.0))/(4.0*pi);

    // intervalo numérico padrão
    double zMin = 1.0e-6;
    double zMax = 1.0 - 1.0e-6;

    // r em GeV^{-1}
    double rMin = 1.0e-6;
    double rMax = 1.0e2;
};

inline double massOf(QuarkFlavor f, const QuarkMasses& masses = QuarkMasses{}) {
    switch (f) {
        case QuarkFlavor::u: return masses.u;
        case QuarkFlavor::d: return masses.d;
        case QuarkFlavor::s: return masses.s;
        case QuarkFlavor::c: return masses.c;
        case QuarkFlavor::b: return masses.b;
        case QuarkFlavor::t: return masses.t;
    }
    throw std::runtime_error("Sabor de quark inválido.");
}

inline ElectroweakCouplings makeCouplings(CurrentType current, QuarkFlavor flavor) {
    if (current == CurrentType::CC) {
        return {
            -1.0,
             1.0,
            (MW*MW*GF/std::sqrt(2.0))/(4.0*pi)
        };
    }

    // NC: Z0 -> q qbar
    const double alphaEW = (MZ*MZ*GF/std::sqrt(2.0))/(4.0*pi);

    switch (flavor) {
        case QuarkFlavor::d:
        case QuarkFlavor::s:
        case QuarkFlavor::b:
            return {
                -0.5 + (2.0/3.0)*sin2ThetaW,  // gV = T3 - 2Q sin^2(thetaW)
                -0.5,                           // gA = T3
                alphaEW
            };

        case QuarkFlavor::u:
        case QuarkFlavor::c:
        case QuarkFlavor::t:
            return {
                0.5 - (4.0/3.0)*sin2ThetaW,   // gV = T3 - 2Q sin^2(thetaW)
                0.5,                            // gA = T3
                alphaEW
            };
    }

    throw std::runtime_error("Sabor de quark inválido em NC.");
}

inline Parameters makeNCParameters(
    QuarkFlavor flavor,
    const QuarkMasses& masses = QuarkMasses{}
) {
    Parameters p;
    p.current = CurrentType::NC;
    p.flavor = flavor;

    p.m = massOf(flavor, masses);
    p.mu = p.m;

    auto c = makeCouplings(CurrentType::NC, flavor);
    p.gV = c.gV;
    p.gA = c.gA;
    p.alphaEW = c.alphaEW;

    return p;
}

inline Parameters makeCCParameters(
    QuarkFlavor quark,
    QuarkFlavor antiquark,
    const QuarkMasses& masses = QuarkMasses{}
) {
    Parameters p;
    p.current = CurrentType::CC;
    p.flavor = quark;

    p.m = massOf(quark, masses);
    p.mu = massOf(antiquark, masses);

    auto c = makeCouplings(CurrentType::CC, quark);
    p.gV = c.gV;
    p.gA = c.gA;
    p.alphaEW = c.alphaEW;

    return p;
}

} // namespace dipole

#endif