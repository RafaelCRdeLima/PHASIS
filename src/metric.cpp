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

double Metric::turning_factor_at_turn(double r_turn) const
{
    // T(r_t) = 2 f(r_t) r_t - f'(r_t) r_t^2,  f' por diferenca central.
    if (!(r_turn > 0.0)) return 0.0;
    const double eps = 1.0e-6*r_turn;
    const double fp  = (f(r_turn + eps) - f(r_turn - eps))/(2.0*eps);
    return 2.0*f(r_turn)*r_turn - fp*r_turn*r_turn;
}

double Metric::turning_factor_slope(double r_turn) const
{
    // T'(r_t) = f(r_t) - (1/2) f''(r_t) r_t^2,  f'' por diferenca central.
    if (!(r_turn > 0.0)) return 1.0;
    const double eps = 1.0e-4*r_turn;
    const double fpp = (f(r_turn + eps) - 2.0*f(r_turn) + f(r_turn - eps))/(eps*eps);
    return f(r_turn) - 0.5*fpp*r_turn*r_turn;
}

double Metric::r_turning(double b, double r_search_max) const
{
    // Ponto de retorno: menor raiz externa de  g(r) = r/sqrt(f(r)) - b.
    //
    // Implementacao generica: procura um bracket descendo a partir de
    // r_search_max e depois bisecta. Metricas com forma fechada devem
    // sobrescrever (Minkowski e Schwarzschild o fazem).
    //
    // O bracketing NAO pode simplesmente comecar num r pequeno: em
    // metricas com horizonte, f < 0 la dentro e sqrt(f) da NaN. Por isso
    // a descida para assim que f deixa de ser positivo.
    //
    // Esta busca ja implementa, de graca, o criterio de captura: em
    // Schwarzschild g tem MINIMO na esfera de fotons, valendo
    // b_crit - b. Para b < b_crit nao ha mudanca de sinal, a descida
    // termina sem bracket, e o metodo lanca -- que e o que
    // Metric::is_captured interpreta como captura.

    if (!(b > 0.0)) return 0.0;

    auto g = [&](double r) { return r/std::sqrt(f(r)) - b; };

    double hi = r_search_max;
    if (!(f(hi) > 0.0)) {
        throw std::domain_error("Metric::r_turning: f <= 0 em r_search_max");
    }
    if (g(hi) < 0.0) {
        throw std::domain_error(
            "Metric::r_turning: sem ponto de retorno ate r_search_max");
    }

    double lo = hi;
    bool achou = false;
    for (int i = 0; i < 400; ++i) {
        const double prox = lo*0.7;
        if (!(f(prox) > 0.0)) break;      // chegou ao horizonte sem cruzar
        lo = prox;
        if (g(lo) < 0.0) { achou = true; break; }
    }

    if (!achou) {
        throw std::domain_error(
            "Metric::r_turning: sem mudanca de sinal acima do horizonte; "
            "o raio e capturado");
    }

    for (int i = 0; i < 300; ++i) {
        const double mid = 0.5*(lo + hi);
        if (g(mid) > 0.0) hi = mid; else lo = mid;
        if (hi - lo <= 1.0e-16*hi) break;
    }
    return 0.5*(lo + hi);
}

// =====================================================================
// Metric: default generico de captura
// =====================================================================
bool Metric::is_captured(double b) const
{
    if (!(b > 0.0)) return false;
    try {
        (void)r_turning(b, b);
        return false;
    } catch (const std::domain_error&) {
        return true;   // sem ponto de retorno = engolido
    }
}

// =====================================================================
// Schwarzschild
// =====================================================================
namespace {
constexpr double kSqrt3 = 1.7320508075688772;
// b_crit / r_s = 3 sqrt(3) / 2
constexpr double kBcritOverRs = 2.598076211353316;
}

Schwarzschild Schwarzschild::from_solar_masses(double M_sol)
{
    return Schwarzschild(2.95325e5*M_sol);
}

double Schwarzschild::b_crit() const
{
    return kBcritOverRs*r_s_;
}

bool Schwarzschild::is_captured(double b) const
{
    // Discriminante da cubica r_t^3 - b^2 r_t + b^2 r_s = 0 e
    //     Delta = b^4 (4 b^2 - 27 r_s^2),
    // positivo <=> b > (3 sqrt3 / 2) r_s. Nao ha nada aproximado aqui:
    // o criterio de captura E o discriminante.
    if (!(b > 0.0)) return true;                 // raio radial cai
    return !(4.0*b*b > 27.0*r_s_*r_s_);
}

double Schwarzschild::r_turning(double b, double) const
{
    if (is_captured(b)) {
        throw std::domain_error(
            "Schwarzschild::r_turning: b <= b_crit, o raio e capturado");
    }

    // arg = -3 sqrt3 r_s / (2b),  em [-1, 0) para b > b_crit.
    double arg = -1.5*kSqrt3*r_s_/b;
    if (arg < -1.0) arg = -1.0;    // so por seguranca de arredondamento
    if (arg >  1.0) arg =  1.0;

    return (2.0*b/kSqrt3)*std::cos(std::acos(arg)/3.0);
}

double Schwarzschild::turning_factor(double r, double r_turn) const
{
    if (!(r_turn > 0.0)) {
        // r_t = 0: T = f(0) r^2 / r, mas f(0) e singular em Schwarzschild.
        throw std::domain_error(
            "Schwarzschild::turning_factor: r_turn = 0 nao tem sentido aqui");
    }

    // Forma exata e fatorada: nenhuma subtracao de numeros proximos.
    return (r + r_turn) - (r_s_/(r_turn*r))*(r*r + r*r_turn + r_turn*r_turn);
}

} // namespace phasis
