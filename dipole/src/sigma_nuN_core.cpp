#include "sigma_nuN_core.hpp"

#include <cmath>

namespace dipole {

double largeXFactor(double x)
{
    if (x >= 1.0) return 0.0;
    if (x <= 0.0) return 0.0;
    return std::pow(1.0 - x, 7.0);
}

// =====================================================
// F4 - numero de nos proporcional a faixa de integracao
//
// O codigo antigo usava NlogQ = 16 FIXO para uma faixa
// ln(1) ate ln(0.999*2*M_N*E), que cresce com a energia. Em 1e6 GeV
// isso da ~5 nos cobrindo o pico em Q^2 ~ M_W^2; em 1e14 GeV, ~2.
// O resultado parava de crescer e ficava ruidoso, e uma perturbacao
// de 1e-7 em E movia sigma em 2%.
//
// Fixando a DENSIDADE de nos por decada, o espacamento em ln Q^2 nao
// muda com a energia e o pico e sempre resolvido igual.
// =====================================================
int nodesForRange(double log_lo, double log_hi,
                         double per_decade, int n_min)
{
    const double decades = (log_hi - log_lo)/std::log(10.0);

    int n = static_cast<int>(std::ceil(per_decade*decades));
    if (n < n_min) n = n_min;
    if (n % 2 != 0) ++n;

    return n;
}

double d2sigma_dxdy_CC(
    double Enu,
    double x,
    double Q2,
    const StructureTable& table,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam
)
{
    const double s = 2.0*MN*Enu;
    const double y = Q2/(x*s);

    if (x <= 0.0 || x >= 1.0) return 0.0;
    if (y <= 0.0 || y >= 1.0) return 0.0;

    // A tabela guarda F_T e F_L crus; o (1-x)^7 (regra de contagem de
    // constituintes, Kutak-Kwiecinski) e aplicado aqui para manter a
    // tabela suave.
    const StructureTL raw = table.at(x, Q2);
    const double lx = largeXFactor(x);

    StructureTL F;
    F.FT = lx*raw.FT;
    F.FL = lx*raw.FL;

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
        0.5*(1.0 + std::pow(1.0-y, 2))*(F.FT + F.FL)
        - 0.5*y*y*F.FL
        + f3_sign*y*(1.0 - 0.5*y)*xF3;

    return prefactor*bracket;
}

double sigmaNuN_CC(
    double Enu,
    const StructureTable& table,
    double nodesPerDecadeQ,
    double nodesPerDecadeX,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam,
    double Q2min,
    int* nodesQout
)
{
    const double s = 2.0*MN*Enu;
    const double Q2max = 0.999*s;

    if (Q2max <= Q2min) return 0.0;

    const double logQmin = std::log(Q2min);
    const double logQmax = std::log(Q2max);

    const int NlogQ = nodesForRange(logQmin, logQmax, nodesPerDecadeQ, 16);
    if (nodesQout) *nodesQout = NlogQ;

    const double xmax = 0.999999;

    // O numero de nos em x tem de ser FIXO ao longo do laco de Q^2.
    //
    // Se for recalculado por no (a faixa e [Q^2/s, xmax], que encolhe
    // conforme Q^2 cresce), ele salta de 2 em 2 ao longo do laco, e cada
    // salto e uma descontinuidade no integrando EXTERNO: o erro de
    // convergencia da integral interna muda de degrau. Deslocar os nos
    // externos por 1e-7 fazia alguns cruzarem esses degraus, e sigma
    // pulava 0.2-0.6% - o que sobrava do ruido diagnosticado na secao 0.
    //
    // Fixando pela faixa mais larga (a de Q^2 = Q2min), o integrando
    // externo volta a ser suave e as faixas estreitas ficam apenas
    // super-resolvidas, o que so custa tempo.
    const int Nlogx = nodesForRange(std::log(Q2min/s), std::log(xmax),
                                    nodesPerDecadeX, 16);

    auto integrandLogQ = [&](double logQ2)
    {
        const double Q2 = std::exp(logQ2);

        const double xmin = Q2/s;

        if (xmin >= xmax) return 0.0;

        const double logxmin = std::log(xmin);
        const double logxmax = std::log(xmax);

        auto integrandLogx = [&](double logx)
        {
            const double x = std::exp(logx);

            return d2sigma_dxdy_CC(
                       Enu, x, Q2, table, useF3, weakSF, beam
                   )/s;
        };

        return Q2*simpson(integrandLogx, logxmin, logxmax, Nlogx);
    };

    return simpson(integrandLogQ, logQmin, logQmax, NlogQ);
}


} // namespace dipole
