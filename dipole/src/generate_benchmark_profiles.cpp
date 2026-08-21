// ================================================================
// Benchmark profiles for DIS optical-depth tests with neutrinos
//
// Output files:
//
//   data/benchmark_01_radiating_sphere.dat
//   data/benchmark_02_collapsar_disk_funnel.dat
//   data/benchmark_03_uhe_attenuation.dat
//
// Columns:
//
//   r_cm dens_gcm3 temp_MeV ye pres_dummy entr_dummy
//
// Notes:
//
// Benchmark 1:
//   Spherical radiating-source atmosphere.
//   Useful as a clean transport/free-streaming/redshift test.
//
// Benchmark 2:
//   Collapsar-like power-law envelope with approximate disk/funnel
//   encoded as an effective radial profile.
//
// Benchmark 3:
//   Ultra-high-energy attenuation target.
//   Dense inner region plus extended envelope, designed to test
//   tau(E) = integral n(r) sigma_DIS(E) dr.
//
// ================================================================

#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>
#include <string>

double clamp(double x, double xmin, double xmax)
{
    if (x < xmin) return xmin;
    if (x > xmax) return xmax;
    return x;
}

void write_row(std::ofstream& out, double r, double rho, double temp, double ye)
{
    const double rho_safe = std::max(rho, 1.0e-99);
    double pres = rho * temp;
    double entr = temp / std::pow(rho_safe, 1.0/3.0);

    out << r << " "
        << rho << " "
        << temp << " "
        << ye << " "
        << pres << " "
        << entr << "\n";
}

int main()
{
    std::filesystem::create_directories("data");

    const int N = 3000;

    const double rmin = 1.0e6;   // 10 km
    const double rmax = 1.0e11;  // 1e6 km

    const double rho_floor = 1.0e-12;
    const double T_floor   = 1.0e-4;

    // ============================================================
    // Benchmark 1: relativistic radiating-sphere-like atmosphere
    // ============================================================
    {
        std::ofstream out("data/benchmark_01_radiating_sphere.dat");

        const double Rsrc = 4.0e5;     // 4 km, source radius
        const double r0   = 1.0e7;     // 100 km
        const double rho0 = 1.0e8;     // g/cm^3
        const double T0   = 8.0;       // MeV
        const double Rcut = 5.0e9;     // cm

        out << "# Benchmark 1: radiating-sphere-like atmosphere\n";
        out << "# rho(r) = rho0 / [1 + (r/r0)^2] * exp[-r/Rcut]\n";
        out << "# T(r)   = T0   / [1 + (r/r0)^0.5]\n";
        out << "# Rsrc_cm = " << Rsrc << "\n";
        out << "# columns: r_cm dens_gcm3 temp_MeV ye pres_dummy entr_dummy\n";

        for (int i = 0; i < N; ++i) {
            double logr = std::log(rmin)
                        + i*(std::log(rmax)-std::log(rmin))/(N-1.0);
            double r = std::exp(logr);

            double rho = rho0 / (1.0 + std::pow(r/r0, 2.0))
                       * std::exp(-r/Rcut);

            rho = std::max(rho, rho_floor);

            double temp = T0 / (1.0 + std::pow(r/r0, 0.5));
            temp = std::max(temp, T_floor);

            double ye = 0.30;

            write_row(out, r, rho, temp, ye);
        }
    }

    // ============================================================
    // Benchmark 2: collapsar-like disk/funnel effective profile
    // ============================================================
    {
        std::ofstream out("data/benchmark_02_collapsar_disk_funnel.dat");

        const double r0       = 1.0e7;   // 100 km
        const double rho0     = 1.0e9;   // g/cm^3 at r0
        const double alpha    = 2.5;
        const double Rcut     = 3.0e9;   // 30000 km
        const double funnel_f = 1.0e-2;  // effective polar dilution
        const double disk_f   = 5.0;     // effective disk enhancement

        const double T0       = 10.0;    // MeV
        const double betaT    = 0.8;

        out << "# Benchmark 2: synthetic collapsar-like disk/funnel effective profile\n";
        out << "# base rho(r) = rho0*(r/r0)^(-alpha)*exp[-(r-r0)/Rcut]\n";
        out << "# effective rho = base rho * sqrt(disk_f * funnel_f)\n";
        out << "# alpha = " << alpha << "\n";
        out << "# disk_f = " << disk_f << ", funnel_f = " << funnel_f << "\n";
        out << "# columns: r_cm dens_gcm3 temp_MeV ye pres_dummy entr_dummy\n";

        for (int i = 0; i < N; ++i) {
            double logr = std::log(rmin)
                        + i*(std::log(rmax)-std::log(rmin))/(N-1.0);
            double r = std::exp(logr);

            double base = rho0
                        * std::pow(r/r0, -alpha)
                        * std::exp(-(r-r0)/Rcut);

            // Effective 1D average between dense disk and low-density polar funnel.
            double geometry_factor = std::sqrt(disk_f * funnel_f);

            double rho = base * geometry_factor;
            rho = std::max(rho, rho_floor);

            double temp = T0
                        * std::pow(r/r0, -betaT)
                        * std::exp(-(r-r0)/(2.0*Rcut));

            temp = std::max(temp, T_floor);

            double ye = 0.10 + 0.40*(std::log(r/rmin)/std::log(rmax/rmin));
            ye = clamp(ye, 0.05, 0.50);

            write_row(out, r, rho, temp, ye);
        }
    }

    // ============================================================
    // Benchmark 3: UHE attenuation / DIS optical-depth profile
    // ============================================================
    {
        std::ofstream out("data/benchmark_03_uhe_attenuation.dat");

        const double rc    = 5.0e6;   // 50 km
        const double rho_c = 1.0e11;  // g/cm^3
        const double n     = 3.0;

        const double rho_env = 1.0e6;   // g/cm^3
        const double r_env   = 1.0e9;   // 10000 km
        const double Rcut    = 1.0e10;  // 100000 km

        const double T_c   = 20.0;    // MeV
        const double rT    = 1.0e7;   // 100 km

        out << "# Benchmark 3: UHE neutrino attenuation / DIS optical-depth profile\n";
        out << "# rho(r) = rho_c/[1+(r/rc)^n] + rho_env*(r/r_env)^(-2)*exp[-r/Rcut]\n";
        out << "# Designed for tau(E) = integral n_b(r)*sigma_DIS(E)*dr tests\n";
        out << "# columns: r_cm dens_gcm3 temp_MeV ye pres_dummy entr_dummy\n";

        for (int i = 0; i < N; ++i) {
            double logr = std::log(rmin)
                        + i*(std::log(rmax)-std::log(rmin))/(N-1.0);
            double r = std::exp(logr);

            double core = rho_c / (1.0 + std::pow(r/rc, n));

            double env = rho_env
                       * std::pow(r/r_env, -2.0)
                       * std::exp(-r/Rcut);

            double rho = core + env;
            rho = std::max(rho, rho_floor);

            double temp = T_c / (1.0 + r/rT);
            temp = std::max(temp, T_floor);

            // More neutron rich in the dense region.
            double ye = 0.08 + 0.32*(std::log(r/rmin)/std::log(rmax/rmin));
            ye = clamp(ye, 0.05, 0.45);

            write_row(out, r, rho, temp, ye);
        }
    }

    std::cout << "Arquivos gerados:\n";
    std::cout << "  data/benchmark_01_radiating_sphere.dat\n";
    std::cout << "  data/benchmark_02_collapsar_disk_funnel.dat\n";
    std::cout << "  data/benchmark_03_uhe_attenuation.dat\n";

    return 0;
}