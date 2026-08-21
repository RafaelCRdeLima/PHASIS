// ================================================================
// Synthetic collapsar-like regularized core-envelope profile
//
// Density profile:
//
//     rho(r) = rho_c / [ 1 + (r/rc)^3 ]
//
// Temperature profile:
//
//     T(r) = T_c / [ 1 + (r/rT) ]
//
// This type of parametrized profile is commonly used in
// neutrino transport, ray tracing, leakage schemes,
// and collapsar-inspired test problems because:
//
//   - it avoids central singularities;
//   - produces a dense proto-compact core;
//   - reproduces a power-law envelope at large radii;
//   - is numerically stable for optical depth integration;
//   - mimics collapsar and post-bounce CCSN atmospheres.
//
// Similar parametrized atmospheres and power-law envelopes
// appear in:
//
// [1] O'Connor, E. (2012)
//     "An Open-Source Neutrino Radiation Hydrodynamics Code for
//      Core-Collapse Supernovae"
//     Astrophysical Journal Supplement, 219, 24.
//
// [2] Sekiguchi, Y. et al. (2012)
//     "Dynamical Mass Ejection from the Merger of Binary Neutron Stars"
//     Progress of Theoretical and Experimental Physics.
//
// [3] Ott, C. D. et al. (2008)
//     "Multiparameter Study of the Gravitational-Wave Signature
//      of Core-Collapse Supernovae"
//     Physical Review D, 78, 084069.
//
// ================================================================


#include <iostream>
#include <fstream>
#include <cmath>

int main()
{
    const int N = 3000;

    // Raio em cm: 10 km até 1e11 cm = 1e6 km
    const double rmin = 1.0e6;
    const double rmax = 1.0e11;

    // Escala de referência: 100 km
    const double r0 = 1.0e7;

    // Perfil tipo collapsar/envelope:
    // denso no centro, queda forte e corte exponencial externo
    const double rho0 = 1.0e9;      // g/cm^3 em r0
    const double alpha = 3.0;       // queda tipo potência
    const double Rcut = 3.0e9;      // cm = 30000 km
    const double rho_floor = 1.0e-12;

    // Temperatura
    const double T0 = 10.0;         // MeV em r0
    const double betaT = 0.8;
    const double T_floor = 1.0e-4;  // MeV

    std::ofstream out("data/collapsar_profile.dat");

    out << "# Collapsar-like analytic profile with exponential cutoff\n";
    out << "# rho(r) = rho0*(r/r0)^(-alpha)*exp[-(r-r0)/Rcut]\n";
    out << "# columns:\n";
    out << "# r_cm dens_gcm3 temp_MeV ye pres_dummy entr_dummy\n";

    for (int i = 0; i < N; ++i) {

        double logr = std::log(rmin)
                    + i*(std::log(rmax)-std::log(rmin))/(N-1.0);

        double r = std::exp(logr);

        double power_part = std::pow(r/r0, -alpha);
        double cutoff_part = std::exp(-(r-r0)/Rcut);

        double rho = rho0 * power_part * cutoff_part;

        if (rho < rho_floor) rho = rho_floor;

        double temp = T0 * std::pow(r/r0, -betaT)
                    * std::exp(-(r-r0)/(2.0*Rcut));

        if (temp < T_floor) temp = T_floor;

        // neutron-rich perto do centro, mais protonizado fora
        double ye = 0.10 + 0.40*(std::log(r/rmin)/std::log(rmax/rmin));
        if (ye > 0.5) ye = 0.5;
        if (ye < 0.05) ye = 0.05;

        // Campos auxiliares apenas para manter formato FLASH-like
        double pres = rho * temp;
        double entr = temp / std::pow(rho, 1.0/3.0);

        out << r << " "
            << rho << " "
            << temp << " "
            << ye << " "
            << pres << " "
            << entr << "\n";
    }

    std::cout << "Arquivo gerado: data/collapsar_profile.dat\n";

    return 0;
}