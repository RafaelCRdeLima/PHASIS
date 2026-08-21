#ifndef DIPOLE_INTEGRALS_HPP
#define DIPOLE_INTEGRALS_HPP

#include <functional>

#include "parameters.hpp"
#include "dipole_models.hpp"

namespace dipole {

// integração simples Simpson
double simpson(
    const std::function<double(double)>& f,
    double a,
    double b,
    int N
);

// integração dupla
double integrate2D(
    const std::function<double(double,double)>& f,
    double ax,
    double bx,
    int Nx,
    double ay,
    double by,
    int Ny
);

// integrandos físicos
double integrandT(
    double r,
    double z,
    double x,
    double Q2,
    const Parameters& wf,
    const GBWParameters& gbw
);

double integrandL(
    double r,
    double z,
    double x,
    double Q2,
    const Parameters& wf,
    const GBWParameters& gbw
);

// funções de estrutura
double FT_GBW(
    double x,
    double Q2,
    const Parameters& wf,
    const GBWParameters& gbw,
    int Nr = 400,
    int Nz = 200
);

double FL_GBW(
    double x,
    double Q2,
    const Parameters& wf,
    const GBWParameters& gbw,
    int Nr = 400,
    int Nz = 200
);

double F2_GBW(
    double x,
    double Q2,
    const Parameters& wf,
    const GBWParameters& gbw,
    int Nr = 400,
    int Nz = 200
);


double integrandT_IIM(
    double r,
    double z,
    double x,
    double Q2,
    const Parameters& wf,
    const IIMParameters& iim
);

double integrandL_IIM(
    double r,
    double z,
    double x,
    double Q2,
    const Parameters& wf,
    const IIMParameters& iim
);

double FT_IIM(
    double x,
    double Q2,
    const Parameters& wf,
    const IIMParameters& iim,
    int Nr = 400,
    int Nz = 200
);

double FL_IIM(
    double x,
    double Q2,
    const Parameters& wf,
    const IIMParameters& iim,
    int Nr = 400,
    int Nz = 200
);

double F2_IIM(
    double x,
    double Q2,
    const Parameters& wf,
    const IIMParameters& iim,
    int Nr = 400,
    int Nz = 200
);

} // namespace dipole

#endif