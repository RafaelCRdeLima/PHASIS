#include "integrals.hpp"

#include "wavefunctions.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include <boost/math/special_functions/bessel.hpp>

namespace dipole {

namespace {

// Nos e pesos de Simpson composto em [a,b] com N intervalos (N par).
struct SimpsonRule {
    std::vector<double> node;
    std::vector<double> weight;
};

SimpsonRule simpsonRule(double a, double b, int N)
{
    if (N < 2) N = 2;
    if (N % 2 != 0) ++N;

    const double h = (b - a)/N;

    SimpsonRule rule;
    rule.node.resize(N + 1);
    rule.weight.resize(N + 1);

    for (int i = 0; i <= N; ++i) {
        rule.node[i] = a + i*h;

        double w;
        if (i == 0 || i == N) w = 1.0;
        else if (i % 2 != 0)  w = 4.0;
        else                  w = 2.0;

        rule.weight[i] = w*h/3.0;
    }

    return rule;
}

} // namespace

double simpson(
    const std::function<double(double)>& f,
    double a,
    double b,
    int N
)
{
    const SimpsonRule rule = simpsonRule(a, b, N);

    double sum = 0.0;
    for (std::size_t i = 0; i < rule.node.size(); ++i) {
        sum += rule.weight[i]*f(rule.node[i]);
    }

    return sum;
}

double integrate2D(
    const std::function<double(double,double)>& f,
    double ax, double bx, int Nx,
    double ay, double by, int Ny
)
{
    auto gx = [&](double x)
    {
        auto gy = [&](double y) { return f(x,y); };
        return simpson(gy, ay, by, Ny);
    };

    return simpson(gx, ax, bx, Nx);
}

// =====================================================
// Nucleo da quadratura
//
//   F_{T,L} = Q^2/(4 pi^2) * int d^2r int_0^1 dz |psi|^2 sigma_dip
//           = Q^2/(4 pi^2) * int dr (2 pi r) int dz |psi|^2 sigma_dip
//
// Em r: substituicao u = ln r, dr = r du  =>  fator r^2.
// Em z: as duas pontas em escala log. Para a ponta de cima,
//       z = 1 - e^t com o MESMO t, o que da o mesmo peso e^t:
//
//   int_0^1 f dz = int_{ln zMin}^{ln 0.5} e^t [ f(e^t) + f(1-e^t) ] dt
//
// Nao se assume simetria z <-> 1-z: ela vale so quando m = mu, e
// falha no canal (c,s).
// =====================================================
StructureTL structureFunctionsTL(
    double x,
    double Q2,
    const Parameters& wf,
    const std::function<double(double,double)>& sigmaDip,
    const QuadratureGrid& grid
)
{
    if (!(Q2 > 0.0)) {
        throw std::runtime_error("Q2 deve ser positivo em structureFunctionsTL.");
    }
    if (!(grid.rMin > 0.0) || !(grid.rMax > grid.rMin)) {
        throw std::runtime_error("Faixa de r invalida em structureFunctionsTL.");
    }
    if (!(grid.zMin > 0.0) || !(grid.zMin < 0.5)) {
        throw std::runtime_error("zMin invalido em structureFunctionsTL.");
    }

    const SimpsonRule ru = simpsonRule(std::log(grid.rMin),
                                       std::log(grid.rMax), grid.Nr);

    const SimpsonRule rz = simpsonRule(std::log(grid.zMin),
                                       std::log(0.5), grid.Nz/2);

    double FT = 0.0;
    double FL = 0.0;

    for (std::size_t i = 0; i < ru.node.size(); ++i) {

        const double r  = std::exp(ru.node[i]);
        const double sd = sigmaDip(r, x);

        if (sd <= 0.0) continue;

        // 2 pi r (do d^2r) * r (do jacobiano dr = r du)
        const double wr = ru.weight[i] * 2.0*pi * r * r * sd;

        for (std::size_t j = 0; j < rz.node.size(); ++j) {

            const double et = std::exp(rz.node[j]);
            const double wz = rz.weight[j] * et;

            // ponta de baixo: z = e^t ; ponta de cima: z = 1 - e^t
            const double zs[2] = { et, 1.0 - et };

            for (int k = 0; k < 2; ++k) {

                const double z = zs[k];

                const double eps2v = epsilon2(z, Q2, wf);
                if (!(eps2v > 0.0)) continue;

                const double arg = std::sqrt(eps2v)*r;

                const double K0 = boost::math::cyl_bessel_k(0, arg);
                const double K1 = boost::math::cyl_bessel_k(1, arg);

                FT += wr * wz * psiT2_pre(z, Q2, wf, eps2v, K0, K1);
                FL += wr * wz * psiL2_pre(z, Q2, wf, eps2v, K0, K1);
            }
        }
    }

    const double norm = Q2/(4.0*pi*pi);

    StructureTL out;
    out.FT = norm*FT;
    out.FL = norm*FL;
    return out;
}

// =====================================================
// Integrandos de diagnostico
// =====================================================

double integrandT(
    double r, double z, double x, double Q2,
    const Parameters& wf, const GBWParameters& gbw
)
{
    return 2.0*pi*r*psiT2(r, z, Q2, wf)*sigmaDipoleGBW(r, x, gbw);
}

double integrandL(
    double r, double z, double x, double Q2,
    const Parameters& wf, const GBWParameters& gbw
)
{
    return 2.0*pi*r*psiL2(r, z, Q2, wf)*sigmaDipoleGBW(r, x, gbw);
}

// =====================================================
// GBW
// =====================================================

StructureTL structureGBW(
    double x, double Q2, const Parameters& wf,
    const GBWParameters& gbw, const QuadratureGrid& grid
)
{
    auto sd = [&](double r, double xx) { return sigmaDipoleGBW(r, xx, gbw); };
    return structureFunctionsTL(x, Q2, wf, sd, grid);
}

double FT_GBW(double x, double Q2, const Parameters& wf,
              const GBWParameters& gbw, int Nr, int Nz)
{
    QuadratureGrid g; g.Nr = Nr; g.Nz = Nz;
    return structureGBW(x, Q2, wf, gbw, g).FT;
}

double FL_GBW(double x, double Q2, const Parameters& wf,
              const GBWParameters& gbw, int Nr, int Nz)
{
    QuadratureGrid g; g.Nr = Nr; g.Nz = Nz;
    return structureGBW(x, Q2, wf, gbw, g).FL;
}

double F2_GBW(double x, double Q2, const Parameters& wf,
              const GBWParameters& gbw, int Nr, int Nz)
{
    QuadratureGrid g; g.Nr = Nr; g.Nz = Nz;
    const StructureTL F = structureGBW(x, Q2, wf, gbw, g);
    return F.FT + F.FL;
}

// =====================================================
// IIM / bCGC
// =====================================================

StructureTL structureIIM(
    double x, double Q2, const Parameters& wf,
    const IIMParameters& iim, const QuadratureGrid& grid
)
{
    auto sd = [&](double r, double xx) { return sigmaDipoleIIM(r, xx, iim); };
    return structureFunctionsTL(x, Q2, wf, sd, grid);
}

double FT_IIM(double x, double Q2, const Parameters& wf,
              const IIMParameters& iim, int Nr, int Nz)
{
    QuadratureGrid g; g.Nr = Nr; g.Nz = Nz;
    return structureIIM(x, Q2, wf, iim, g).FT;
}

double FL_IIM(double x, double Q2, const Parameters& wf,
              const IIMParameters& iim, int Nr, int Nz)
{
    QuadratureGrid g; g.Nr = Nr; g.Nz = Nz;
    return structureIIM(x, Q2, wf, iim, g).FL;
}

double F2_IIM(double x, double Q2, const Parameters& wf,
              const IIMParameters& iim, int Nr, int Nz)
{
    QuadratureGrid g; g.Nr = Nr; g.Nz = Nz;
    const StructureTL F = structureIIM(x, Q2, wf, iim, g);
    return F.FT + F.FL;
}

} // namespace dipole
