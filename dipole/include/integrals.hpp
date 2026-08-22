#ifndef DIPOLE_INTEGRALS_HPP
#define DIPOLE_INTEGRALS_HPP

#include <functional>

#include "parameters.hpp"
#include "dipole_models.hpp"

namespace dipole {

// =====================================================
// Grade de quadratura
//
// O integrando de F_{T,L} tem duas estruturas que uma grade
// uniforme nao resolve:
//
//   * em r, o pico fica em r ~ 1/eps, com eps ~ sqrt(z(1-z)Q^2).
//     Em Q^2 = M_W^2 isso vale ~0.025 GeV^-1, enquanto a grade
//     uniforme antiga [1e-6, 1e2] com Nr=50 tinha passo 2 GeV^-1:
//     o pico inteiro caia entre o primeiro e o segundo no.
//
//   * em z, ha quase-singularidades integraveis em z->0 e z->1
//     (configuracoes aligned-jet, onde eps -> massa do quark e o
//     dipolo fica grande). Uma grade uniforme da peso espurio ao
//     no de borda.
//
// Os dois erros tem sinais opostos e se cancelam parcialmente:
// corrigir so um PIORA o resultado. Por isso log em r e as duas
// pontas de z resolvidas entram juntos.
//
// Ver CAMPANHA_CORRECAO.md, secao 0.
// =====================================================

struct QuadratureGrid {
    // r em GeV^{-1}, integrado em ln r
    double rMin = 1.0e-8;
    double rMax = 1.0e3;

    // z integrado em ln z a partir das DUAS pontas
    double zMin = 1.0e-11;

    int Nr = 200;   // intervalos de Simpson em ln r
    int Nz = 200;   // intervalos totais em z (metade por ponta)
};

struct StructureTL {
    double FT = 0.0;
    double FL = 0.0;
};

// integração simples Simpson
double simpson(
    const std::function<double(double)>& f,
    double a,
    double b,
    int N
);

// integração dupla (grade uniforme; mantida para uso didatico)
double integrate2D(
    const std::function<double(double,double)>& f,
    double ax,
    double bx,
    int Nx,
    double ay,
    double by,
    int Ny
);

// =====================================================
// Nucleo: F_T e F_L numa unica varredura da grade.
//
// sigmaDip(r, x) e a seção de choque dipolo-proton do modelo.
// K0 e K1 sao calculados uma vez por ponto e usados nas duas.
// =====================================================
StructureTL structureFunctionsTL(
    double x,
    double Q2,
    const Parameters& wf,
    const std::function<double(double,double)>& sigmaDip,
    const QuadratureGrid& grid
);

// integrandos físicos (diagnostico; usados por main.cpp --mode integrand)
double integrandT(
    double r, double z, double x, double Q2,
    const Parameters& wf, const GBWParameters& gbw
);

double integrandL(
    double r, double z, double x, double Q2,
    const Parameters& wf, const GBWParameters& gbw
);

// =====================================================
// GBW
// =====================================================
StructureTL structureGBW(
    double x, double Q2, const Parameters& wf,
    const GBWParameters& gbw, const QuadratureGrid& grid = QuadratureGrid{}
);

double FT_GBW(double x, double Q2, const Parameters& wf,
              const GBWParameters& gbw, int Nr = 200, int Nz = 200);

double FL_GBW(double x, double Q2, const Parameters& wf,
              const GBWParameters& gbw, int Nr = 200, int Nz = 200);

double F2_GBW(double x, double Q2, const Parameters& wf,
              const GBWParameters& gbw, int Nr = 200, int Nz = 200);

// =====================================================
// IIM / bCGC
// =====================================================
StructureTL structureIIM(
    double x, double Q2, const Parameters& wf,
    const IIMParameters& iim, const QuadratureGrid& grid = QuadratureGrid{}
);

double FT_IIM(double x, double Q2, const Parameters& wf,
              const IIMParameters& iim, int Nr = 200, int Nz = 200);

double FL_IIM(double x, double Q2, const Parameters& wf,
              const IIMParameters& iim, int Nr = 200, int Nz = 200);

double F2_IIM(double x, double Q2, const Parameters& wf,
              const IIMParameters& iim, int Nr = 200, int Nz = 200);

} // namespace dipole

#endif
