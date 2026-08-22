#include "phasis/ode.hpp"

#include <algorithm>
#include <cmath>

namespace phasis {

namespace {

// Tableau de Dormand-Prince 5(4). FSAL: k7 do passo e o k1 do proximo.
constexpr double c2 = 1.0/5.0,  c3 = 3.0/10.0, c4 = 4.0/5.0,
                 c5 = 8.0/9.0,  c6 = 1.0,      c7 = 1.0;

constexpr double a21 =  1.0/5.0;
constexpr double a31 =  3.0/40.0,        a32 =  9.0/40.0;
constexpr double a41 = 44.0/45.0,        a42 = -56.0/15.0,      a43 = 32.0/9.0;
constexpr double a51 = 19372.0/6561.0,   a52 = -25360.0/2187.0,
                 a53 = 64448.0/6561.0,   a54 = -212.0/729.0;
constexpr double a61 = 9017.0/3168.0,    a62 = -355.0/33.0,
                 a63 = 46732.0/5247.0,   a64 = 49.0/176.0,      a65 = -5103.0/18656.0;
constexpr double a71 = 35.0/384.0,       a73 = 500.0/1113.0,    a74 = 125.0/192.0,
                 a75 = -2187.0/6784.0,   a76 = 11.0/84.0;

// solucao de 4a ordem, para a estimativa de erro
constexpr double b1s = 5179.0/57600.0,   b3s = 7571.0/16695.0,  b4s = 393.0/640.0,
                 b5s = -92097.0/339200.0, b6s = 187.0/2100.0,   b7s = 1.0/40.0;

} // namespace

template <std::size_t N>
OdeStats dopri54(const std::function<void(double, const OdeState<N>&, OdeState<N>&)>& f,
                 double s0, double s1, OdeState<N>& y, const OdeOpts& opts)
{
    OdeStats st;

    const double span = s1 - s0;
    if (!(span > 0.0)) return st;

    OdeState<N> k1, k2, k3, k4, k5, k6, k7, tmp, y5, y4;

    f(s0, y, k1);
    ++st.n_evals;

    double h = (opts.h_init > 0.0) ? opts.h_init : span*1.0e-3;
    const double h_min = opts.h_min_rel*span;
    const double h_max = opts.h_max_rel*span;
    h = std::min(h, h_max);
    double s = s0;
    bool fsal = true;

    while (s < s1) {
        if (st.steps >= opts.max_steps) {
            st.ok = false; st.hit_max_steps = true; break;
        }

        h = std::min(h, s1 - s);
        if (h < h_min) { st.ok = false; st.hit_min_step = true; break; }

        if (!fsal) { f(s, y, k1); ++st.n_evals; }

        for (std::size_t i = 0; i < N; ++i) tmp[i] = y[i] + h*(a21*k1[i]);
        f(s + c2*h, tmp, k2);

        for (std::size_t i = 0; i < N; ++i) tmp[i] = y[i] + h*(a31*k1[i] + a32*k2[i]);
        f(s + c3*h, tmp, k3);

        for (std::size_t i = 0; i < N; ++i)
            tmp[i] = y[i] + h*(a41*k1[i] + a42*k2[i] + a43*k3[i]);
        f(s + c4*h, tmp, k4);

        for (std::size_t i = 0; i < N; ++i)
            tmp[i] = y[i] + h*(a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
        f(s + c5*h, tmp, k5);

        for (std::size_t i = 0; i < N; ++i)
            tmp[i] = y[i] + h*(a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
        f(s + c6*h, tmp, k6);

        for (std::size_t i = 0; i < N; ++i)
            y5[i] = y[i] + h*(a71*k1[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
        f(s + c7*h, y5, k7);

        st.n_evals += 6;

        for (std::size_t i = 0; i < N; ++i)
            y4[i] = y[i] + h*(b1s*k1[i] + b3s*k3[i] + b4s*k4[i]
                            + b5s*k5[i] + b6s*k6[i] + b7s*k7[i]);

        // erro relativo por componente, norma do maximo
        double err = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            const double sc = opts.abs_tol
                            + opts.rel_tol*std::max(std::fabs(y[i]), std::fabs(y5[i]));
            if (sc > 0.0) err = std::max(err, std::fabs(y5[i] - y4[i])/sc);
        }

        if (err <= 1.0 || h <= h_min) {
            s += h;
            y = y5;
            k1 = k7;          // FSAL
            fsal = true;
            ++st.steps;
        } else {
            ++st.rejected;
            fsal = false;
        }

        // controlador padrao, expoente 1/5, com limites de seguranca
        const double fator = (err > 0.0)
            ? std::min(5.0, std::max(0.2, 0.9*std::pow(1.0/err, 0.2)))
            : 5.0;
        h = std::min(h*fator, h_max);
    }

    st.h_final = h;
    return st;
}

// ---------------------------------------------------------------------
// Versao de dimensao dinamica. Mesmo tableau, mesmo controlador; a unica
// diferenca e o piso por componente e a alocacao dos k_i.
OdeStats dopri54_dyn(
    const std::function<void(double, const std::vector<double>&, std::vector<double>&)>& f,
    double s0, double s1, std::vector<double>& y, const OdeOpts& opts)
{
    OdeStats st;

    const double span = s1 - s0;
    if (!(span > 0.0)) return st;

    const std::size_t N = y.size();
    const bool por_comp = (opts.abs_tol_per_component.size() == N);

    std::vector<double> k1(N), k2(N), k3(N), k4(N), k5(N), k6(N), k7(N),
                        tmp(N), y5(N), y4(N);

    f(s0, y, k1);
    ++st.n_evals;

    double h = (opts.h_init > 0.0) ? opts.h_init : span*1.0e-3;
    const double h_min = opts.h_min_rel*span;
    const double h_max = opts.h_max_rel*span;
    h = std::min(h, h_max);
    double s = s0;
    bool fsal = true;

    while (s < s1) {
        if (st.steps >= opts.max_steps) { st.ok = false; st.hit_max_steps = true; break; }

        h = std::min(h, s1 - s);
        if (h < h_min) { st.ok = false; st.hit_min_step = true; break; }

        if (!fsal) { f(s, y, k1); ++st.n_evals; }

        for (std::size_t i = 0; i < N; ++i) tmp[i] = y[i] + h*(a21*k1[i]);
        f(s + c2*h, tmp, k2);
        for (std::size_t i = 0; i < N; ++i) tmp[i] = y[i] + h*(a31*k1[i] + a32*k2[i]);
        f(s + c3*h, tmp, k3);
        for (std::size_t i = 0; i < N; ++i)
            tmp[i] = y[i] + h*(a41*k1[i] + a42*k2[i] + a43*k3[i]);
        f(s + c4*h, tmp, k4);
        for (std::size_t i = 0; i < N; ++i)
            tmp[i] = y[i] + h*(a51*k1[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
        f(s + c5*h, tmp, k5);
        for (std::size_t i = 0; i < N; ++i)
            tmp[i] = y[i] + h*(a61*k1[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
        f(s + c6*h, tmp, k6);
        for (std::size_t i = 0; i < N; ++i)
            y5[i] = y[i] + h*(a71*k1[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
        f(s + c7*h, y5, k7);
        st.n_evals += 6;

        for (std::size_t i = 0; i < N; ++i)
            y4[i] = y[i] + h*(b1s*k1[i] + b3s*k3[i] + b4s*k4[i]
                            + b5s*k5[i] + b6s*k6[i] + b7s*k7[i]);

        double err = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            const double piso = por_comp ? opts.abs_tol_per_component[i] : opts.abs_tol;
            const double sc = piso
                            + opts.rel_tol*std::max(std::fabs(y[i]), std::fabs(y5[i]));
            if (sc > 0.0) err = std::max(err, std::fabs(y5[i] - y4[i])/sc);
        }

        if (err <= 1.0 || h <= h_min) {
            s += h;
            y = y5;
            k1 = k7;
            fsal = true;
            ++st.steps;
        } else {
            ++st.rejected;
            fsal = false;
        }

        const double fator = (err > 0.0)
            ? std::min(5.0, std::max(0.2, 0.9*std::pow(1.0/err, 0.2)))
            : 5.0;
        h = std::min(h*fator, h_max);
    }

    st.h_final = h;
    return st;
}

template OdeStats dopri54<6>(
    const std::function<void(double, const OdeState<6>&, OdeState<6>&)>&,
    double, double, OdeState<6>&, const OdeOpts&);

} // namespace phasis
