// ================================================================
// F2 - teste de aceitacao da quadratura
//
// 1) Valor de referencia: F_2 do canal (u,d) em x=1e-5, Q^2=M_W^2.
//    Convergido = 18.0 com F_L/F_2 = 0.099.
//    A grade uniforme antiga (Nr=50, Nz=30) dava 1018 com F_L/F_2 = 4e-5.
//
// 2) Convergencia: dobrar Nr e Nz deve mudar F_2 em < 1%.
//
// 3) Normalizacao: a soma dos dois canais, com massas -> 0 e dividida
//    por alphaEW, deve reproduzir Kutak-Kwiecinski EPJ C29 (2003) 521,
//    eqs. (12) e (13).
// ================================================================

#include <cmath>
#include <cstdio>
#include <vector>

#include "parameters.hpp"
#include "dipole_models.hpp"
#include "integrals.hpp"
#include "wavefunctions.hpp"

using namespace dipole;

namespace {

const double F2_REF     = 18.0;
const double FL_F2_REF  = 0.099;

int falhas = 0;

void checa(const char* nome, double got, double esperado, double tol_rel)
{
    const double rel = std::fabs(got - esperado)/std::fabs(esperado);
    const bool ok = rel <= tol_rel;
    if (!ok) ++falhas;
    std::printf("  %-52s %12.5g  (ref %.5g, desvio %6.2f%%)  %s\n",
                nome, got, esperado, 100.0*rel, ok ? "OK" : "<-- FALHA");
}

StructureTL canalUD(double x, double Q2, int Nr, int Nz)
{
    Parameters wf = makeCCParameters(QuarkFlavor::u, QuarkFlavor::d);
    GBWParameters gbw;
    QuadratureGrid g; g.Nr = Nr; g.Nz = Nz;
    StructureTL F = structureGBW(x, Q2, wf, gbw, g);
    F.FT /= wf.alphaEW;
    F.FL /= wf.alphaEW;
    return F;
}

} // namespace

int main()
{
    const double x  = 1.0e-5;
    const double Q2 = MW*MW;

    std::printf("\n=== F2: teste de aceitacao da quadratura ===\n\n");
    std::printf("1) Valor de referencia  (canal u,d;  x = 1e-5,  Q2 = M_W^2)\n\n");

    const StructureTL F = canalUD(x, Q2, 200, 200);
    checa("F_2", F.FT + F.FL, F2_REF, 0.05);
    checa("F_L/F_2", F.FL/(F.FT + F.FL), FL_F2_REF, 0.10);

    std::printf("\n2) Convergencia (dobrar a grade muda < 1%%)\n\n");

    const int grades[][2] = { {100,100}, {200,200}, {400,400}, {800,800} };
    std::vector<double> vals;
    for (const auto& g : grades) {
        const StructureTL Fg = canalUD(x, Q2, g[0], g[1]);
        const double f2 = Fg.FT + Fg.FL;
        vals.push_back(f2);
        std::printf("  Nr = %-4d Nz = %-4d   F_2 = %12.6f   F_L/F_2 = %.5f\n",
                    g[0], g[1], f2, Fg.FL/f2);
    }
    std::printf("\n");
    for (std::size_t i = 1; i < vals.size(); ++i) {
        const double rel = std::fabs(vals[i] - vals[i-1])/vals[i-1];
        const bool ok = rel <= 0.01;
        if (!ok) ++falhas;
        std::printf("  variacao %zu->%zu: %7.4f%%   %s\n",
                    i-1, i, 100.0*rel, ok ? "OK" : "<-- FALHA");
    }

    std::printf("\n3) Normalizacao vs Kutak-Kwiecinski eqs. (12) e (13)\n");
    std::printf("   (soma dos 2 canais, massas -> 0, dividido por alphaEW)\n\n");

    {
        QuarkMasses leves;
        leves.u = leves.d = leves.s = leves.c = 1.0e-9;

        const double rr = 1.0, zz = 0.5, QQ = 10.0;
        const double Qb2 = zz*(1.0 - zz)*QQ;

        // KK eq (12) e (13)
        const double kkT = 6.0/(pi*pi)*(zz*zz + (1-zz)*(1-zz))*Qb2
                         * std::pow(std::cyl_bessel_k(1, std::sqrt(Qb2)*rr), 2);
        const double kkL = 24.0/(pi*pi)*zz*zz*(1-zz)*(1-zz)*QQ
                         * std::pow(std::cyl_bessel_k(0, std::sqrt(Qb2)*rr), 2);

        double codT = 0.0, codL = 0.0;
        const QuarkFlavor pares[2][2] = {
            { QuarkFlavor::u, QuarkFlavor::d },
            { QuarkFlavor::c, QuarkFlavor::s }
        };
        for (const auto& par : pares) {
            Parameters wf = makeCCParameters(par[0], par[1], leves);
            const double e2 = epsilon2(zz, QQ, wf);
            const double ar = std::sqrt(e2)*rr;
            codT += psiT2_pre(zz, QQ, wf, e2, std::cyl_bessel_k(0, ar),
                              std::cyl_bessel_k(1, ar)) / wf.alphaEW;
            codL += psiL2_pre(zz, QQ, wf, e2, std::cyl_bessel_k(0, ar),
                              std::cyl_bessel_k(1, ar)) / wf.alphaEW;
        }

        checa("|psi_T|^2  vs KK eq (12)", codT, kkT, 1.0e-6);
        checa("|psi_L|^2  vs KK eq (13)", codL, kkL, 1.0e-6);
    }

    std::printf("\n  RESULTADO: %s\n\n", falhas == 0 ? "OK" : "FALHA");
    return falhas == 0 ? 0 : 1;
}
