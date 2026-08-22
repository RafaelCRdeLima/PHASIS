#ifndef DIPOLE_WAVEFUNCTIONS_HPP
#define DIPOLE_WAVEFUNCTIONS_HPP

#include <cmath>

#include "parameters.hpp"

namespace dipole {

// K_0 e K_1 num ponto, com guarda no argumento grande.
//
// std::cyl_bessel_k (C++17) LANCA excecao para argumento muito grande
// -- no libstdc++, em algum ponto entre 1e6 e 1e9 -- enquanto o boost
// simplesmente devolvia zero por underflow. Como
// K_nu(x) ~ sqrt(pi/2x) e^{-x}, em x = 700 ja vale menos que 1e-305:
// multiplicado por qualquer fator do integrando (o maior e eps^2, que
// chega a ~5e13 na faixa tabelada) o produto e exatamente zero em dupla
// precisao. Devolver zero acima do limiar reproduz o boost e evita a
// excecao.
//
// Isso importa porque a tabela cobre Q^2 ate ~2e14 GeV^2 e r ate
// 1e3 GeV^-1, o que da eps*r ate ~1e9.
inline void besselK01(double arg, double& K0, double& K1)
{
    constexpr double LIMIAR = 700.0;

    if (!std::isfinite(arg) || arg > LIMIAR) {
        K0 = 0.0;
        K1 = 0.0;
        return;
    }

    K0 = std::cyl_bessel_k(0, arg);
    K1 = std::cyl_bessel_k(1, arg);
}

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
