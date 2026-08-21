#include <iostream>
#include <fstream>
#include <cmath>

int main()
{
    const int N = 1000;

    // Raio em cm
    const double rmin = 0.0;
    const double rmax = 2.0e8; // 2000 km

    // Perfil analítico tipo proto-NS / CCSN
    const double rho_c = 3.0e14; // g/cm^3
    const double Rrho  = 1.0e6;  // cm
    const double n     = 8.0;

    const double T_c   = 30.0;   // MeV
    const double RT    = 2.0e6;  // cm

    std::ofstream out("data/flash_like_profile.dat");

    out << "# FLASH-like analytic profile\n";
    out << "# columns:\n";
    out << "# r_cm dens_gcm3 temp_MeV ye pres_dummy entr_dummy\n";

    for (int i = 0; i < N; ++i) {

        double r = rmin + (rmax - rmin)*i/(N - 1.0);

        double rho = rho_c / (1.0 + std::pow(r/Rrho, n));

        double temp = T_c * std::exp(-r/RT) + 0.1;

        double ye = 0.05 + 0.45*(r/rmax);
        if (ye > 0.5) ye = 0.5;

        // Campos extras apenas para imitar tabela hidrodinâmica
        double pres = rho * temp;   // dummy
        double entr = temp / std::pow(rho, 1.0/3.0); // dummy

        out << r << " "
            << rho << " "
            << temp << " "
            << ye << " "
            << pres << " "
            << entr << "\n";
    }

    std::cout << "Arquivo gerado: data/flash_like_profile.dat\n";

    return 0;
}