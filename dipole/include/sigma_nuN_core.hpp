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

double d2sigma_dxdy_CC(
    double Enu, double x, double Q2,
    const StructureTable& table,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam
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

} // namespace dipole

#endif
