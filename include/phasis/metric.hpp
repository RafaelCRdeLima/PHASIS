#ifndef PHASIS_METRIC_HPP
#define PHASIS_METRIC_HPP

#include <string>

namespace phasis {

// =====================================================================
// Metrica estatica e esfericamente simetrica:
//
//     ds^2 = -f(r) dt^2 + h(r) dr^2 + r^2 dOmega^2
//
// Convencao de assinatura (-,+,+,+).
// =====================================================================
struct Metric {
    virtual ~Metric() = default;

    virtual double f(double r) const = 0;   // -g_tt
    virtual double h(double r) const = 0;   //  g_rr

    virtual std::string name() const = 0;

    // Menor raiz de  b = r / sqrt(f(r)),  o ponto de retorno do raio de
    // parametro de impacto b.
    //
    // A implementacao default e numerica (bisecao). Metricas com forma
    // fechada devem sobrescrever -- e o que garante os 1e-10 dos testes
    // T1/T2, onde qualquer erro de raiz apareceria direto no tau.
    //
    // Lanca std::domain_error se nao houver ponto de retorno na faixa
    // buscada. A semantica de CAPTURA (buraco negro) entra na Fase 2.
    virtual double r_turning(double b, double r_search_max) const;

    // T(r) = [ f(r_t) r^2 - f(r) r_t^2 ] / (r - r_t),   r_t = ponto de retorno.
    //
    // Por que isto existe. O integrando tem o fator 1/sqrt(w) com
    //
    //     w(r) = 1 - f(r) b^2 / r^2 ,     b^2 = r_t^2 / f(r_t)
    //
    // e w(r_t) = 0. Calcular w por essa formula perto de r_t e subtrair
    // dois numeros quase iguais: em r - r_t ~ 1e-16 r nao sobra digito
    // algum. Reescrevendo,
    //
    //     w(r) = (r - r_t) T(r) / ( f(r_t) r^2 )
    //
    // e, com a substituicao r = r_t + s^2 (dr = 2s ds), o s CANCELA
    // analiticamente:
    //
    //     dl/dr * 2s = 2 r sqrt( h(r) f(r_t) / T(r) )
    //
    // Sem singularidade e sem cancelamento. T e a "parte suave" de w, e
    // cada metrica deve fornece-la em forma fechada quando puder.
    virtual double turning_factor(double r, double r_turn) const;

    // T no proprio ponto de retorno. A forma geral de turning_factor da
    // 0/0 em r = r_t, entao o limite entra separado:
    //
    //     T(r_t) = 2 f(r_t) r_t - f'(r_t) r_t^2
    //
    // Esta quantidade e o diagnostico central do regime quase-critico:
    // ela mede a distancia a raiz DUPLA de w. Em Schwarzschild vale
    // 2 r_t - 3 r_s e zera na esfera de fotons.
    virtual double turning_factor_at_turn(double r_turn) const;

    // T'(r_t) = f(r_t) - (1/2) f''(r_t) r_t^2.  Da 1 tanto em Minkowski
    // quanto em Schwarzschild; serve para dimensionar o pico do
    // integrando, s_* = sqrt( T(r_t) / T'(r_t) ).
    virtual double turning_factor_slope(double r_turn) const;

    // Raio da esfera de fotons: o r onde g(r) = r/sqrt(f(r)) tem MINIMO.
    //
    //     g'(r) = (1/2) r^(1/2) (r-r_s)^(-3/2) (2r - 3 r_s)
    //
    // Zero em r = 1.5 r_s para Schwarzschild. Metricas sem minimo (g
    // monotonica, como Minkowski) devolvem 0.
    //
    // Nao e curiosidade: e o que decide o SENTIDO do criterio de escape
    // de um raio emitido. Abaixo da esfera de fotons g DECRESCE, entao um
    // raio emitido para FORA encontra g(r) = b acima de r_emit, vira, e
    // cai. A condicao de escape inverte.
    virtual double r_photon_sphere() const { return 0.0; }

    // Existe ponto de retorno para este b?
    //
    // ATENCAO: isto NAO e "o raio e capturado" em geral. Para um raio
    // vindo do infinito as duas coisas coincidem; para um raio EMITIDO
    // nao. Ver phasis::classify em emission.hpp.
    //
    // Default generico: nao ha raiz externa, ou seja, r_turning falha.
    // Metricas que sabem o criterio em forma fechada devem sobrescrever
    // -- em Schwarzschild o criterio E o discriminante da cubica, o que
    // torna o teste exato em vez de "o solver nao convergiu".
    virtual bool is_captured(double b) const;

    // Energia local medida por um observador estatico em r, para um raio
    // com energia conservada E_inf = f * dt/dlambda.
    //
    //     E_loc = -p_mu u^mu = E_inf / sqrt(f(r))
    //
    // Em Minkowski isto e a identidade; ja fica escrito assim para a
    // Fase 2 nao precisar tocar aqui.
    double E_local(double E_inf, double r) const;
};

// ---------------------------------------------------------------------
struct Minkowski final : Metric {
    double f(double) const override { return 1.0; }
    double h(double) const override { return 1.0; }
    std::string name() const override { return "Minkowski"; }

    // f = 1  =>  b = r  =>  r_turning = b, exato.
    double r_turning(double b, double) const override { return b; }

    // f = 1:  T = (r^2 - r_t^2)/(r - r_t) = r + r_t.  Exato, sem
    // subtracao de termos proximos.
    double turning_factor(double r, double r_turn) const override {
        return r + r_turn;
    }
    double turning_factor_at_turn(double r_turn) const override { return 2.0*r_turn; }
    double turning_factor_slope(double)          const override { return 1.0; }
};

// ---------------------------------------------------------------------
// Schwarzschild:  f(r) = 1/h(r) = 1 - r_s/r.
//
// r_s = 2GM/c^2 = 2.95325e5 cm * (M/M_sol).
// ---------------------------------------------------------------------
// Nao e `final`: a suite de testes deriva dela para forcar a bisecao
// generica e validar as duas rotas de r_turning uma contra a outra.
class Schwarzschild : public Metric {
public:
    explicit Schwarzschild(double r_s_cm) : r_s_(r_s_cm) {}

    static Schwarzschild from_solar_masses(double M_sol);

    double f(double r) const override { return 1.0 - r_s_/r; }
    double h(double r) const override { return 1.0/(1.0 - r_s_/r); }
    std::string name() const override { return "Schwarzschild"; }

    // Forma fechada. b = r_t/sqrt(1 - r_s/r_t) e a cubica
    //
    //     r_t^3 - b^2 r_t + b^2 r_s = 0
    //
    // com discriminante b^4 (4 b^2 - 27 r_s^2). A raiz fisica (a maior,
    // k = 0 na solucao trigonometrica) e
    //
    //     r_t = (2b/sqrt3) cos[ (1/3) acos( -3 sqrt3 r_s / (2b) ) ]
    //
    // Exata, sem iteracao. Remove o root-finder como fonte de erro na
    // metrica de referencia. Limites: r_s -> 0 da r_t = b;
    // b = b_crit da r_t = 1.5 r_s.
    double r_turning(double b, double r_search_max) const override;

    // T(r) = (r + r_t) - (r_s/(r_t r)) (r^2 + r r_t + r_t^2)
    //
    // vem de  f(r_t) r^2 - f(r) r_t^2 = (r^2 - r_t^2) - r_s (r^3 - r_t^3)/(r_t r),
    // fatorando (r - r_t) dos dois termos. Exata.
    //
    // No ponto de retorno:  T(r_t) = 2 r_t - 3 r_s.
    // Isso zera exatamente na esfera de fotons, r_ph = 1.5 r_s -- nao e
    // coincidencia: la r_t vira raiz DUPLA de w, e e por isso que existe
    // b_crit. Ver turning_factor_at_turn().
    double turning_factor(double r, double r_turn) const override;

    // T(r_t) = 2 r_t - 3 r_s, sem passar pela forma geral (onde r = r_t
    // faria 0/0).
    double turning_factor_at_turn(double r_turn) const override {
        return 2.0*r_turn - 3.0*r_s_;
    }

    // T'(r_t) = f(r_t) - (1/2) f''(r_t) r_t^2
    //         = (1 - r_s/r_t) - (1/2)(-2 r_s/r_t^3) r_t^2 = 1.  Exato.
    double turning_factor_slope(double) const override { return 1.0; }

    // Captura <=> discriminante da cubica <= 0 <=> b <= (3 sqrt3 / 2) r_s.
    bool is_captured(double b) const override;

    double r_s() const { return r_s_; }
    double b_crit() const;              // (3 sqrt(3)/2) r_s
    double r_photon() const { return 1.5*r_s_; }
    double r_photon_sphere() const override { return 1.5*r_s_; }

private:
    double r_s_;
};

} // namespace phasis

#endif
