#ifndef PHASIS_EMISSION_HPP
#define PHASIS_EMISSION_HPP

#include <string>

#include "phasis/metric.hpp"
#include "phasis/ray.hpp"

namespace phasis {

// =====================================================================
// Emissor ESTATICO em r_emit.
//
// No referencial ortonormal estatico local, com psi medido da direcao
// radial PARA FORA:
//
//     L      = r_emit E_loc sin(psi)
//     E_inf  = sqrt(f(r_emit)) E_loc
//     b      = L/E_inf = r_emit sin(psi) / sqrt(f(r_emit))
//
// Auto-verificacao: em psi = pi/2 isto da b = r/sqrt(f) = g(r), que e
// exatamente a relacao do ponto de retorno -- um raio emitido
// tangencialmente tem r_t = r_emit.
//
// Emissao isotropica no referencial local significa cos(psi) uniforme em
// [-1,1], NAO psi uniforme.
//
// Emissor ORBITANDO exige aberracao adicional entre o referencial do
// fluido e o estatico, e nao esta implementado. O efeito nao e pequeno:
// beta orbital = sqrt(M/r)/sqrt(1-2M/r) vale EXATAMENTE 1/2 na ISCO
// (r = 6M), com gamma = 2/sqrt(3), e o fator Doppler gamma(1 +- beta)
// varre 1/sqrt(3) a sqrt(3) -- fator 3 em energia, contra o teto
// sqrt(3) do redshift gravitacional. Adia-lo e sequenciamento, nao a
// afirmacao de que e pequeno.
// =====================================================================

double b_from_emission(double r_emit, double psi, const Metric& metric);

// ---------------------------------------------------------------------
// Topologia do raio emitido.
//
// g(r) = r/sqrt(f(r)) tem MINIMO b_crit na esfera de fotons, entao o
// criterio de escape depende de que lado dela r_emit esta -- sao QUATRO
// casos, nao tres:
//
//   r_emit > r_ph,  psi < pi/2  ->  OutboundOnly   (sempre escapa)
//   r_emit > r_ph,  psi > pi/2  ->  b > b_crit ? Turning : Captured
//   r_emit < r_ph,  psi < pi/2  ->  b < b_crit ? OutboundOnly : Captured
//                                   (INVERTIDA: o raio sobe, vira acima
//                                    de r_emit, e cai)
//   r_emit < r_ph,  psi > pi/2  ->  Captured (sempre)
//
// Em r_emit = r_ph os dois ramos concordam (b <= b_crit sempre, com
// igualdade so em psi = pi/2), entao nao ha descontinuidade.
// ---------------------------------------------------------------------
enum class Topology { OutboundOnly, Turning, Captured };

const char* to_string(Topology t);

Topology classify(double r_emit, double psi, const Metric& metric,
                  double* b_out = nullptr);

// Angulo do cone de escape:
//     sin(psi_c) = b_crit sqrt(f(r)) / r
// e fracao capturada com cos(psi) uniforme:
//     f_cap = (1 - cos(psi_c))/2
//
// Dois valores exatos:  r = 1.5 r_s -> psi_c = pi/2, f_cap = 1/2
//                       r = 3.0 r_s -> psi_c = pi/4, f_cap = (2-sqrt2)/4
double escape_cone_angle(double r_emit, const Metric& metric);
double captured_fraction(double r_emit, const Metric& metric);

// Constroi um raio a partir do angulo de emissao.
Ray emit_ray(double r_emit, double psi, double E_inf_GeV,
             const Metric& metric,
             double inclination_rad = 0.0, double psi_turn_rad = 0.0);

} // namespace phasis

#endif
