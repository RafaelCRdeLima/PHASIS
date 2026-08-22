#include "phasis/metric.hpp"

#include <cmath>
#include <stdexcept>

namespace phasis {

double Metric::E_local(double E_inf, double r) const
{
    const double fr = f(r);
    if (!(fr > 0.0)) {
        throw std::domain_error("Metric::E_local: f(r) <= 0 (dentro do horizonte?)");
    }
    return E_inf/std::sqrt(fr);
}

double Metric::turning_factor(double r, double r_turn) const
{
    // Default generico. Longe de r_turn a divisao direta e estavel;
    // perto dela cai para o limite
    //
    //     T(r_t) = lim = 2 f(r_t) r_t - f'(r_t) r_t^2
    //
    // com f' por diferenca central. Metricas com forma fechada devem
    // sobrescrever -- Minkowski ja o faz, e Schwarzschild fara na Fase 2
    // (la vale  T = (r + r_t) - (r_s/(r_t r))(r^2 + r r_t + r_t^2) ).

    const double f_t  = f(r_turn);
    const double d    = r - r_turn;
    const double escala = (r_turn > 0.0) ? r_turn : r;

    if (std::fabs(d) > 1.0e-6*escala) {
        return (f_t*r*r - f(r)*r_turn*r_turn)/d;
    }

    if (!(r_turn > 0.0)) {
        // r_t = 0: T(r) = f(0) r^2 / r = f(0) r
        return f_t*r;
    }

    const double eps = 1.0e-6*r_turn;
    const double fp  = (f(r_turn + eps) - f(r_turn - eps))/(2.0*eps);
    return 2.0*f_t*r_turn - fp*r_turn*r_turn;
}

double Metric::r_turning(double b, double r_search_max) const
{
    // Ponto de retorno: menor raiz de  g(r) = r/sqrt(f(r)) - b = 0.
    //
    // Implementacao default, generica e deliberadamente simples: bisecao
    // no intervalo (0, r_search_max]. Para f <= 1 vale r/sqrt(f) >= r,
    // logo a raiz satisfaz r_turn <= b, e passar r_search_max = b basta.
    //
    // NAO trata captura: em Schwarzschild g(r) tem MINIMO na esfera de
    // fotons, e para b abaixo do critico nao ha raiz alguma. Detectar
    // isso, e usar Brent em vez de bisecao, e trabalho da Fase 2.
    // Metricas com forma fechada devem sobrescrever este metodo.

    if (!(b > 0.0)) return 0.0;

    auto g = [&](double r) { return r/std::sqrt(f(r)) - b; };

    double lo = r_search_max*1.0e-12;
    double hi = r_search_max;

    if (g(hi) < 0.0) {
        throw std::domain_error(
            "Metric::r_turning: sem ponto de retorno ate r_search_max");
    }
    if (g(lo) > 0.0) {
        throw std::domain_error(
            "Metric::r_turning: sem mudanca de sinal; possivel captura "
            "(tratamento e da Fase 2)");
    }

    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5*(lo + hi);
        if (g(mid) > 0.0) hi = mid; else lo = mid;
        if (hi - lo <= 1.0e-15*hi) break;
    }
    return 0.5*(lo + hi);
}

} // namespace phasis
