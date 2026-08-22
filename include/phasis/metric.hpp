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
};

} // namespace phasis

#endif
