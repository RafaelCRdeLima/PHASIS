#ifndef DIPOLE_SIGMA_NUN_CORE_HPP
#define DIPOLE_SIGMA_NUN_CORE_HPP

#include "structure_table.hpp"
#include "weak_structure_functions.hpp"

namespace dipole {

// 1 GeV^{-2} = 0.389379 mb = 0.389379e-27 cm^2
constexpr double GeVminus2_to_cm2 = 0.389379e-27;
constexpr double MN = 0.938272;   // massa do nucleon em GeV

// Regra de contagem de constituintes, n_s = 4 espectadores.
// Kutak-Kwiecinski EPJ C29 (2003) 521.
double largeXFactor(double x);

// F4: numero de nos proporcional a faixa, para manter a DENSIDADE
// por decada constante. Ver CAMPANHA_CORRECAO.md secao 0, Sentinela 2.
int nodesForRange(double log_lo, double log_hi, double per_decade, int n_min);

// =====================================================================
// FORMULA MESTRA, UMA SO PARA AS DUAS CORRENTES
//
//   d2sigma/dx dy = (G_F^2 M_N E / pi) ( M_V^2/(Q^2+M_V^2) )^2
//                   [ (1+(1-y)^2)/2 F_2 - y^2/2 F_L + s3 y(1-y/2) xF_3 ]
//
// A UNICA diferenca entre CC e NC aqui e M_V (M_W ou M_Z) e de onde vem
// xF_3. Toda a diferenca de acoplamento ja esta DENTRO de F_2 e F_L: a
// tabela divide por alphaEW do proprio canal (convencao DIS), entao o
// que sobra em F e a soma efetiva de cargas -- 4 para o CC (2 canais x
// (g_V^2+g_A^2) = 2), e Soma_q (g_V^2+g_A^2) = 1.3127 para o NC sobre
// u, d, s, c.
//
// Escrever as duas correntes como UMA funcao nao e economia de linhas:
// e a garantia de que elas nao possam divergir. A razao sigma_NC/sigma_CC
// so significa alguma coisa se o resto for identico, e aqui e o mesmo
// codigo.
//
// s3 = +1 para nu, -1 para nubar.
// =====================================================================

// Parametrizada por (x, Q^2); y = Q^2/(x s) sai daqui.
double d2sigma_dxdQ2_param(
    double Enu, double x, double Q2,
    CurrentType current,
    const StructureTable& table,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam
);

// Parametrizada por (x, y); Q^2 = x y s sai daqui. E a forma natural
// para dsigma/dy, e evita o ida-e-volta Q^2 -> y -> Q^2, que perderia
// digitos no no de borda -- justo onde y = 1 e o integrando importa.
double d2sigma_dxdy_param(
    double Enu, double x, double y,
    CurrentType current,
    const StructureTable& table,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam
);

// Compatibilidade: o caminho CC de sempre.
double d2sigma_dxdy_CC(
    double Enu, double x, double Q2,
    const StructureTable& table,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam
);

// =====================================================================
// Secao de choque total, por integracao em (Q^2, x).
// =====================================================================
double sigmaNuN(
    double Enu,
    CurrentType current,
    const StructureTable& table,
    double nodesPerDecadeQ,
    double nodesPerDecadeX,
    bool useF3 = false,
    const weak::WeakStructureFunctions* weakSF = nullptr,
    weak::BeamType beam = weak::BeamType::Neutrino,
    double Q2min = 1.0,
    int* nodesQout = nullptr
);

double sigmaNuN_CC(
    double Enu,
    const StructureTable& table,
    double nodesPerDecadeQ,
    double nodesPerDecadeX,
    bool useF3 = false,
    const weak::WeakStructureFunctions* weakSF = nullptr,
    weak::BeamType beam = weak::BeamType::Neutrino,
    double Q2min = 1.0,
    int* nodesQout = nullptr
);

// =====================================================================
// dsigma/dy a y FIXO, integrando so em x.
//
//   dsigma/dy = INT dx  d2sigma/dx dy,   com Q^2 = x y s
//
// O corte inferior em x vem de Q^2 >= Q2min, ou seja x >= Q2min/(y s):
// abaixo disso o modelo de dipolo nao e valido e a tabela de F nem
// cobre. Isso torna dsigma/dy IDENTICAMENTE NULA para
// y < Q2min/(s * x_max) -- um corte cinematico duro na ponta de baixo
// em y, que quem monta a grade em log y precisa conhecer.
//
// Ordenar a mesma integral dupla de duas maneiras e um teste de
// consistencia forte:  INT dy (dsigma/dy)  tem de dar sigma.
// =====================================================================
double dsigma_dy(
    double Enu,
    double y,
    CurrentType current,
    const StructureTable& table,
    double nodesPerDecadeX,
    bool useF3 = false,
    const weak::WeakStructureFunctions* weakSF = nullptr,
    weak::BeamType beam = weak::BeamType::Neutrino,
    double Q2min = 1.0
);

// Menor y com dsigma/dy nao nula: Q2min/(s * x_max).
double yMinKinematic(double Enu, double Q2min, double xmax = 0.999999);

} // namespace dipole

#endif
