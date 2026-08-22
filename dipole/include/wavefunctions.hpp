#ifndef DIPOLE_WAVEFUNCTIONS_HPP
#define DIPOLE_WAVEFUNCTIONS_HPP

#include "parameters.hpp"

namespace dipole {

double epsilon2(double z, double Q2, const Parameters& p);

double psiT2(double r, double z, double Q2, const Parameters& p);

double psiL2(double r, double z, double Q2, const Parameters& p);

// Versoes que recebem eps2, K0 e K1 ja calculados.
//
// K0 e K1 sao o custo dominante do integrando, e psiT2/psiL2 precisam
// exatamente dos mesmos dois valores no mesmo ponto (r,z). Calcular uma
// vez e passar para as duas permite obter F_T e F_L numa unica varredura
// da grade, em vez de duas integrais 2D independentes.
double psiT2_pre(double z, double Q2, const Parameters& p,
                 double eps2v, double K0, double K1);

double psiL2_pre(double z, double Q2, const Parameters& p,
                 double eps2v, double K0, double K1);

} // namespace dipole

#endif
