#ifndef PHASIS_TRACE_HPP
#define PHASIS_TRACE_HPP

#include "phasis/cross_section.hpp"
#include "phasis/density.hpp"
#include "phasis/integrate.hpp"
#include "phasis/metric.hpp"
#include "phasis/ray.hpp"

namespace phasis {

// =====================================================================
// Profundidade optica ao longo de UMA geodesica nula.
//
// FUNCAO PURA: sem estado global, sem I/O, sem memoria entre chamadas.
// Tudo o que ela toca vem por referencia const. E seguro chamar em
// paralelo -- um `#pragma omp parallel for` por cima e legitimo sem
// nenhuma mudanca aqui -- desde que as implementacoes de Metric,
// DensityProfile e CrossSection passadas tambem sejam const-thread-safe
// (as desta fase sao).
//
//     tau = integral de  N_A * rho(r,theta) * sigma(E_loc(r)) * dl
//
// com
//
//     dl/dr    = sqrt(h) / sqrt(1 - f b^2 / r^2)
//     E_loc(r) = E_inf / sqrt(f(r))
//
// sigma fica DENTRO da integral. Na Fase 1, com f = 1, E_loc e constante
// e a integral fatoraria; nao fatoramos de proposito, porque na Fase 2 o
// redshift faz sigma variar ao longo do raio e nada aqui deve mudar.
// =====================================================================
// Derivadas do raio em relacao a s, com r = r_turn + s^2.
//
//     dl/ds   = 2 r     sqrt( h(r)      f(r_t) / T(r) )
//     dpsi/ds = (2 b/r) sqrt( h(r) f(r) f(r_t) / T(r) )
//
// Exposto para que o teste T13 verifique a identidade
//
//     dpsi/dl = b sqrt(f(r)) / r^2
//
// contra ESTE codigo, o mesmo que a EDO usa -- e nao contra uma
// re-derivacao, que poderia repetir o mesmo erro de fator.
struct RayDerivatives {
    double dl_ds   = 0.0;
    double dpsi_ds = 0.0;
};

RayDerivatives ray_derivatives(const Metric& metric,
                               double r, double r_turn, double b);

Result trace_ray(const Ray& ray,
                 const Metric& metric,
                 const DensityProfile& profile,
                 const CrossSection& xsec,
                 const IntegratorOpts& opts = IntegratorOpts{},
                 // Congela sigma em E_inf, ignorando o redshift no seu
                 // argumento. Mesma assinatura e mesmo significado que
                 // em transport_cascade.
                 //
                 // A GEOMETRIA e a MESMA: o caminho, o ponto de retorno
                 // e a densidade nao mudam, so o argumento de sigma.
                 // Comparar contra Minkowski trocaria o caminho tambem,
                 // e a diferenca medida deixaria de ser so do redshift.
                 //
                 // Com isso,
                 //     R = (tau_full - tau_frozen)/tau_frozen
                 //       = <f^(-alpha/2)> - 1
                 // e uma media ponderada por n(r) dl ao longo do raio --
                 // INVARIANTE sob reescala global de sigma, que e o que
                 // T39 verifica e o que legitima comparar dipolo com
                 // colinear apesar da supressao global do primeiro.
                 bool freeze_redshift = false);

} // namespace phasis

#endif
