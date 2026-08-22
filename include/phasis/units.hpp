#ifndef PHASIS_UNITS_HPP
#define PHASIS_UNITS_HPP

// =====================================================================
// Sistema de unidades do PHASIS.
//
// REGRA UNICA, sem excecao: dentro do codigo tudo esta em
//
//     comprimento .... cm
//     densidade ...... g/cm^3
//     secao de choque. cm^2
//     energia ........ GeV
//     angulo ......... rad
//
// Constantes de conversao existem para ENTRADA e SAIDA (arquivos de
// configuracao, CSV). Nunca para uso interno.
//
// Toda grandeza derivada deve ser checavel dimensionalmente a olho:
//   tau = [1/mol] * [g/cm^3] * [cm^2] * [cm] * [mol/g]  ->  adimensional
// =====================================================================

namespace phasis {
namespace units {

// ---- comprimento ----------------------------------------------------
inline constexpr double cm  = 1.0;
inline constexpr double m   = 100.0 * cm;
inline constexpr double km  = 1.0e5 * cm;

// ---- constantes fisicas ---------------------------------------------

// Numero de Avogadro, por mol.
inline constexpr double N_A = 6.02214076e23;

// Massa molar do nucleon, em g/mol. Usar 1.0 significa tratar
// N_A * rho como densidade de nucleons por cm^3, que e a aproximacao
// padrao em transporte de neutrinos (o desvio e da ordem da energia
// de ligacao nuclear por nucleon, ~1%).
inline constexpr double molar_mass_nucleon_g = 1.0;

// n_nucleon [1/cm^3] = N_A / M_mol * rho [g/cm^3]
inline constexpr double nucleons_per_gram = N_A / molar_mass_nucleon_g;

// ---- constantes para a Fase 2 (metrica) -----------------------------
// r_s = 2GM/c^2. Para o Sol, 2.9532500770e5 cm.
inline constexpr double r_s_per_solar_mass_cm = 2.9532500770e5;

// ---- referencia para testes -----------------------------------------
inline constexpr double earth_radius_cm       = 6371.0 * km;
inline constexpr double earth_mean_density    = 5.51;   // g/cm^3

} // namespace units
} // namespace phasis

#endif
