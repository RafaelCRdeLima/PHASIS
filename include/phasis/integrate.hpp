#ifndef PHASIS_INTEGRATE_HPP
#define PHASIS_INTEGRATE_HPP

#include <functional>

namespace phasis {

struct IntegratorOpts {
    double rel_tol   = 1.0e-8;
    double abs_tol   = 0.0;    // 0 = so tolerancia relativa
    int    max_depth = 50;
};

struct QuadResult {
    double value    = 0.0;
    long   n_evals  = 0;
    int    max_depth_used = 0;
    bool   depth_exhausted = false;   // atingiu max_depth em algum ramo
};

// Simpson adaptativo, sem dependencia externa.
//
// Bisseccao recursiva com o criterio de Richardson |S2 - S1|/15 < tol.
// Retorna tambem o numero de avaliacoes e se algum ramo esgotou a
// profundidade -- o segundo importa: um resultado que esgotou a
// profundidade nao satisfez a tolerancia, e o chamador precisa saber.
QuadResult adaptive_simpson(const std::function<double(double)>& fn,
                            double a, double b,
                            const IntegratorOpts& opts);

} // namespace phasis

#endif
