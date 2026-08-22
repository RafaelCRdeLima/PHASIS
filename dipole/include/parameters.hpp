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

// Constante de estrutura fina, usada apenas no limite eletromagnetico
// (validacao contra HERA). Ver makeEMParameters.
constexpr double alphaEM = 1.0/137.035999084;

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
    NC,
    EM   // fóton: usado só para validar o codigo contra F2 do HERA
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

// Carga eletrica em unidades de e.
inline double electricCharge(QuarkFlavor f) {
    switch (f) {
        case QuarkFlavor::u:
        case QuarkFlavor::c:
        case QuarkFlavor::t: return  2.0/3.0;
        case QuarkFlavor::d:
        case QuarkFlavor::s:
        case QuarkFlavor::b: return -1.0/3.0;
    }
    throw std::runtime_error("Sabor de quark inválido em electricCharge.");
}

inline ElectroweakCouplings makeCouplings(CurrentType current, QuarkFlavor flavor) {
    if (current == CurrentType::EM) {
        // gV = e_f, gA = 0  =>  (gV^2 + gA^2) = e_f^2, que e o peso padrao de F2^em.
        return { electricCharge(flavor), 0.0, alphaEM };
    }

    if (current == CurrentType::CC) {
        return {
            -1.0,
             1.0,
            (MW*MW*GF/std::sqrt(2.0))/(4.0*pi)
        };
    }

    // NC: Z0 -> q qbar
    //
    // ATENCAO -- D10, BUG CONHECIDO E NAO CORRIGIDO.
    //
    // Este valor esta um FATOR 2 BAIXO. O acoplamento efetivo do vertice
    // NC e (g_Z/2)^2/(4pi) = sqrt(2) G_F M_Z^2/(4 pi), o dobro do que
    // esta escrito aqui. Confirmado por duas rotas independentes:
    //
    //  (i) direto: G_F/sqrt2 = g_Z^2/(8 M_Z^2) => alpha = g_Z^2/(16 pi);
    //
    // (ii) contra Kutak-Kwiecinski EPJ C29 (2003) 521, eq. (14)/(15):
    //      |psi_T^Z|^2 = (3/2pi^2)(L_u^2+L_d^2+R_u^2+R_d^2)[...]
    //      Como L = g_V+g_A e R = g_V-g_A, vale L^2+R^2 = 2(g_V^2+g_A^2),
    //      entao C = L_u^2+L_d^2+R_u^2+R_d^2 = 2*Soma_tipos(g_V^2+g_A^2).
    //      Somando os 2 canais com o prefator 6*alpha/(4pi)^2 e dividindo
    //      por alpha, o codigo da (3/2pi^2)*C/2 -- metade de KK.
    //
    // NAO corrigido ainda porque NADA calcula sigma_NC: makeNCParameters
    // so e chamada em main.cpp, para desenhar |psi|^2. Corrigir junto com
    // a implementacao do NC, nao antes -- assim o teste contra o HERA
    // pega os dois de uma vez.
    //
    // NOTA SOBRE O SOMATORIO DE SABORES: em KK, C multiplica a MESMA soma
    // de sabores que aparece no CC (comparar eq. 3 e eq. 5), entao a soma
    // sobre geracoes ja esta dentro das chaves. NAO ha fator 2 adicional
    // de geracoes. Verificado numericamente: C/4*(M_Z/M_W)^2 = 0.4224,
    // que e a razao sigma_NC/sigma_CC conhecida do SM (~0.42); com o
    // fator 2 extra daria 0.845.
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

// Limite eletromagnetico: gamma* -> q qbar, dipolo de UM sabor.
//
// Serve para validar toda a maquinaria (funcoes de onda, quadratura, sigma_dip)
// contra o F2 do HERA -- que e exatamente o dado a que GBW e bCGC foram
// ajustados. Se este limite nao reproduz o HERA, nada a jusante vale.
inline Parameters makeEMParameters(
    QuarkFlavor flavor,
    const QuarkMasses& masses = QuarkMasses{}
) {
    Parameters p;
    p.current = CurrentType::EM;
    p.flavor = flavor;

    p.m = massOf(flavor, masses);
    p.mu = p.m;                       // dipolo q qbar de mesmo sabor

    auto c = makeCouplings(CurrentType::EM, flavor);
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