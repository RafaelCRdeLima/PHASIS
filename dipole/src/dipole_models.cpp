#include "dipole_models.hpp"

#include <cmath>
#include <stdexcept>

namespace dipole {

double mbToGeVminus2(double sigma_mb)
{
    return sigma_mb / 0.389379;
}

// =====================================================
// GBW — parâmetros usados na tese do Alex
// =====================================================

double R0sq_GBW(double x, const GBWParameters& p)
{
    if (x <= 0.0) {
        throw std::runtime_error("x deve ser positivo em R0sq_GBW.");
    }

    return (1.0 / p.Q0sq) * std::pow(x / p.x0, p.lambda);
}

double sigmaDipoleGBW(double r, double x, const GBWParameters& p)
{
    const double sigma0 = mbToGeVminus2(p.sigma0_mb);
    const double R0sq = R0sq_GBW(x, p);

    return sigma0 * (1.0 - std::exp(-r*r/(4.0*R0sq)));
}

// =====================================================
// bCGC / IIM com dependência em parâmetro de impacto
// Versão usada na tese do Alex: forma IIM + Q_s(x,b)
// =====================================================

double Qs_bCGC(double x, double b, const IIMParameters& p)
{
    if (x <= 0.0) {
        throw std::runtime_error("x deve ser positivo em Qs_bCGC.");
    }

    // Q_s^2(x,b) = (x0/x)^lambda * [exp(-b^2/(2 B_CGC))]^(1/gamma_s)
    const double Qs2 =
        std::pow(p.x0 / x, p.lambda)
        * std::pow(
            std::exp(-b*b/(2.0*p.BCGC)),
            1.0/p.gamma_s
        );

    return std::sqrt(Qs2);
}

double amplitude_bCGC(double r, double x, double b, const IIMParameters& p)
{
    if (r <= 0.0) return 0.0;

    const double Y = std::log(1.0/x);
    const double Qs = Qs_bCGC(x, b, p);
    const double tau = r*Qs;

    if (tau <= 0.0) return 0.0;

    if (tau <= 2.0) {
        const double gamma_eff =
            p.gamma_s
            + std::log(2.0/tau)/(9.9*p.lambda*Y);

        double N =
            p.N0 * std::pow(tau/2.0, 2.0*gamma_eff);

        if (N < 0.0) N = 0.0;
        if (N > 1.0) N = 1.0;

        return N;
    }

    // Coeficientes A e B por continuidade em tau = 2.
    const double oneMinusN0 = 1.0 - p.N0;

    const double A =
        - (p.N0*p.N0*p.gamma_s*p.gamma_s)
        / (
            oneMinusN0*oneMinusN0
            * std::log(oneMinusN0)
        );

    const double B =
        0.5 * std::pow(
            oneMinusN0,
            -oneMinusN0/(p.N0*p.gamma_s)
        );

    double N =
        1.0 - std::exp(
            -A * std::pow(std::log(B*tau), 2)
        );

    if (N < 0.0) N = 0.0;
    if (N > 1.0) N = 1.0;

    return N;
}

double sigmaDipoleIIM(double r, double x, const IIMParameters& p)
{
    // sigma_dip(x,r) = 2 ∫ d²b N(x,r,b)
    //                 = 4π ∫ b db N(x,r,b)

    // F7: NAO ha corte em x aqui.
    //
    // O corte duro `x >= 1e-2 -> 0` que existia antes era um desvio da
    // referencia, e era o que produzia a assimetria GBW/IIM de 275 em
    // 1e3 GeV contra 3.3 em 1e14: o GBW extrapolava ate x -> 1 enquanto
    // o IIM zerava.
    //
    // Kutak-Kwiecinski tratam o regime de x grande multiplicando as
    // funcoes de estrutura por (1-x)^(2 n_s - 1) com n_s = 4 espectadores
    // -- a regra de contagem de constituintes, aplicada em
    // sigma_nuN_core.cpp:largeXFactor. A mesma politica vale para os dois
    // modelos. Ver CAMPANHA_CORRECAO.md secao 4, Q4.
    if (x <= 0.0 || x >= 1.0) {
        return 0.0;
    }

    int Nb = p.Nb;
    if (Nb % 2 != 0) ++Nb;

    const double bMin = 0.0;
    const double bMax = p.bMax;
    const double h = (bMax - bMin)/Nb;

    auto integrand = [&](double b)
    {
        return b * amplitude_bCGC(r, x, b, p);
    };

    double sum = integrand(bMin) + integrand(bMax);

    for (int i = 1; i < Nb; ++i) {
        const double b = bMin + i*h;

        if (i % 2 == 0) {
            sum += 2.0*integrand(b);
        } else {
            sum += 4.0*integrand(b);
        }
    }

    const double integral = sum*h/3.0;

    return 4.0*M_PI*integral;
}

} // namespace dipole