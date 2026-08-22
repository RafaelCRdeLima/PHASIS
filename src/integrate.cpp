#include "phasis/integrate.hpp"

#include <cmath>

namespace phasis {

namespace {

struct Ctx {
    const std::function<double(double)>* fn;
    long   n_evals = 0;
    int    max_depth_used = 0;
    bool   exhausted = false;
    int    max_depth;
    double abs_tol;
};

double eval(Ctx& c, double x)
{
    ++c.n_evals;
    return (*c.fn)(x);
}

// Passo recursivo do Simpson adaptativo.
//
// S1 = Simpson em [a,b];  S2 = Simpson nas duas metades.
// Criterio de Richardson: o erro de S2 e ~ |S2 - S1| / 15.
double step(Ctx& c, double a, double b,
            double fa, double fm, double fb, double S1,
            double tol, int depth)
{
    const double m  = 0.5*(a + b);
    const double lm = 0.5*(a + m);
    const double rm = 0.5*(m + b);

    const double flm = eval(c, lm);
    const double frm = eval(c, rm);

    const double h = b - a;
    const double Sl = (h/12.0)*(fa  + 4.0*flm + fm);
    const double Sr = (h/12.0)*(fm  + 4.0*frm + fb);
    const double S2 = Sl + Sr;

    if (depth > c.max_depth_used) c.max_depth_used = depth;

    const double err = std::fabs(S2 - S1);

    if (depth >= c.max_depth) {
        if (err > 15.0*tol) c.exhausted = true;
        return S2 + (S2 - S1)/15.0;
    }

    if (err <= 15.0*tol) {
        // extrapolacao de Richardson: ganha duas ordens de graca
        return S2 + (S2 - S1)/15.0;
    }

    return step(c, a, m, fa, flm, fm, Sl, 0.5*tol, depth + 1)
         + step(c, m, b, fm, frm, fb, Sr, 0.5*tol, depth + 1);
}

} // namespace

QuadResult adaptive_simpson(const std::function<double(double)>& fn,
                            double a, double b,
                            const IntegratorOpts& opts)
{
    QuadResult out;

    if (!(b > a)) return out;   // intervalo vazio ou invertido: zero

    Ctx c;
    c.fn        = &fn;
    c.max_depth = opts.max_depth;
    c.abs_tol   = opts.abs_tol;

    const double m  = 0.5*(a + b);
    const double fa = eval(c, a);
    const double fm = eval(c, m);
    const double fb = eval(c, b);

    const double S1 = ((b - a)/6.0)*(fa + 4.0*fm + fb);

    // Escala para a tolerancia relativa. Se a estimativa inicial for
    // nula (integrando identicamente zero na amostra grosseira), cai
    // para a tolerancia absoluta -- e se essa tambem for zero, usa um
    // piso, para nao exigir precisao infinita de um zero.
    double tol = opts.abs_tol;
    const double scale = std::fabs(S1);
    if (scale > 0.0) tol = std::max(tol, opts.rel_tol*scale);
    if (tol <= 0.0)  tol = 1.0e-300;

    out.value = step(c, a, b, fa, fm, fb, S1, tol, 0);

    out.n_evals         = c.n_evals;
    out.max_depth_used  = c.max_depth_used;
    out.depth_exhausted = c.exhausted;
    return out;
}

} // namespace phasis
