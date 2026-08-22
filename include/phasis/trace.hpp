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
Result trace_ray(const Ray& ray,
                 const Metric& metric,
                 const DensityProfile& profile,
                 const CrossSection& xsec,
                 const IntegratorOpts& opts = IntegratorOpts{});

} // namespace phasis

#endif
