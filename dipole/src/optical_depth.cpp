#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <sstream>

// =====================================================
// Data structures
// =====================================================

struct ProfilePoint {
    double r_cm;
    double rho_gcm3;
    double temp_MeV;
    double ye;
    double pres_dummy;
    double entr_dummy;
};

struct SigmaPoint {
    double Enu_GeV;
    double sigma_GeVminus2;
    double sigma_cm2;
};

// =====================================================
// Readers
// =====================================================

std::vector<ProfilePoint> readProfile(const std::string& filename)
{
    std::ifstream in(filename);

    if (!in) {
        throw std::runtime_error("Could not open profile file: " + filename);
    }

    std::vector<ProfilePoint> profile;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        ProfilePoint p;

        std::istringstream iss(line);

        if (iss >> p.r_cm
                >> p.rho_gcm3
                >> p.temp_MeV
                >> p.ye
                >> p.pres_dummy
                >> p.entr_dummy) {
            profile.push_back(p);
        }
    }

    if (profile.size() < 2) {
        throw std::runtime_error("Profile file has too few points: " + filename);
    }

    return profile;
}

std::vector<SigmaPoint> readSigma(const std::string& filename)
{
    std::ifstream in(filename);

    if (!in) {
        throw std::runtime_error("Could not open sigma table file: " + filename);
    }

    std::vector<SigmaPoint> sigmas;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        SigmaPoint s;

        std::istringstream iss(line);

        // A tabela tem TRES colunas: Enu, sigma[GeV^-2], sigma[cm^2].
        // Ler so duas pegava a coluna em GeV^-2 e a tratava como cm^2,
        // inflando tau por 1/3.894e-28 = 2.6e27.
        if (iss >> s.Enu_GeV >> s.sigma_GeVminus2 >> s.sigma_cm2) {
            sigmas.push_back(s);
        }
    }

    if (sigmas.empty()) {
        throw std::runtime_error("Sigma table is empty: " + filename);
    }

    return sigmas;
}

// =====================================================
// Optical depth
// =====================================================

double computeOpticalDepthTotal(
    const std::vector<ProfilePoint>& profile,
    double sigma_cm2
)
{
    constexpr double m_b = 1.67262192369e-24; // g

    double tau = 0.0;

    for (std::size_t i = 0; i + 1 < profile.size(); ++i) {
        const double r1 = profile[i].r_cm;
        const double r2 = profile[i + 1].r_cm;

        const double dr = std::abs(r2 - r1);

        const double nb1 = profile[i].rho_gcm3 / m_b;
        const double nb2 = profile[i + 1].rho_gcm3 / m_b;

        const double kappa1 = nb1 * sigma_cm2;
        const double kappa2 = nb2 * sigma_cm2;

        tau += 0.5 * (kappa1 + kappa2) * dr;
    }

    return tau;
}

// =====================================================
// Main
// =====================================================

int main(int argc, char* argv[])
{
    std::string profile_file = "data/flash_like_profile.dat";
    std::string sigma_file   = "data/sigma_nuN_CC_GBW.dat";
    std::string output_file  = "data/optical_depth_CC_GBW.dat";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--profile" && i + 1 < argc) {
            profile_file = argv[++i];
        }
        else if (arg == "--sigma-table" && i + 1 < argc) {
            sigma_file = argv[++i];
        }
        else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        }
        else if (arg == "--help" || arg == "-h") {
            std::cout
                << "Usage:\n"
                << "  optical_depth "
                << "--profile data/profile.dat "
                << "--sigma-table data/sigma.dat "
                << "--output data/optical_depth.dat\n";
            return 0;
        }
        else {
            std::cerr << "Unknown or incomplete argument: " << arg << "\n";
            return 1;
        }
    }

    try {
        auto profile = readProfile(profile_file);
        auto sigmas  = readSigma(sigma_file);

        std::ofstream out(output_file);

        if (!out) {
            throw std::runtime_error("Could not open output file: " + output_file);
        }

        out << "# Total optical depth for UHE neutrinos\n";
        out << "# profile_file = " << profile_file << "\n";
        out << "# sigma_file   = " << sigma_file << "\n";
        out << "#\n";
        out << "# Columns:\n";
        out << "# Enu_GeV sigma_cm2 tau_total P_escape\n";

        for (const auto& s : sigmas) {
            const double tau = computeOpticalDepthTotal(profile, s.sigma_cm2);
            const double P_escape = std::exp(-tau);

            out << s.Enu_GeV << " "
                << s.sigma_cm2 << " "
                << tau << " "
                << P_escape << "\n";
        }

        std::cout << "Optical-depth file generated: " << output_file << "\n";
        std::cout << "Profile used: " << profile_file << "\n";
        std::cout << "Sigma table used: " << sigma_file << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error in optical_depth.cpp: " << e.what() << "\n";
        return 1;
    }

    return 0;
}