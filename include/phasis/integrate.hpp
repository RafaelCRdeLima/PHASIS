#ifndef PHASIS_INTEGRATE_HPP
#define PHASIS_INTEGRATE_HPP

#include <functional>

namespace phasis {

struct IntegratorOpts {
    double rel_tol   = 1.0e-8;
    double abs_tol   = 0.0;    // 0 = so tolerancia relativa
    int    max_depth = 50;

    // Orcamento GLOBAL de avaliacoes. Nao e redundante com max_depth:
    // a recursao e binaria, entao limitar a profundidade em 50 permite
    // ate 2^50 folhas. Se a tolerancia pedida cair abaixo do epsilon de
    // maquina (|S2-S1| nao consegue mais encolher), o criterio nunca e
    // satisfeito e a recursao se desdobra ate travar. Este teto e o que
    // transforma isso em resultado degradado e sinalizado, em vez de
    // pendurar o processo.
    long   max_evals = 2000000;

    // Calcular tambem o angulo de deflexao (integral extra, de r_t ate
    // o infinito, independente do perfil de materia).
    bool   want_deflection = false;

    // Limiar de eps = T(r_t)/r_t abaixo do qual o raio e tratado como
    // quase-critico. Ver trace.cpp.
    double near_critical_eps = 1.0e-6;
};

struct QuadResult {
    double value    = 0.0;
    long   n_evals  = 0;
    int    max_depth_used = 0;
    bool   depth_exhausted = false;   // atingiu max_depth em algum ramo
    bool   budget_exhausted = false;  // estourou max_evals
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
