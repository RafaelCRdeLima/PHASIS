#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdexcept>

constexpr double NA = 6.02214076e23;

struct ProfilePoint {
    double r_cm, rho, temp, ye, pres, entr;
};

struct SigmaPoint {
    double Enu_GeV, sigma_GeVminus2, sigma_cm2;
};

std::vector<ProfilePoint> readProfile(const std::string& filename)
{
    std::ifstream in(filename);
    std::vector<ProfilePoint> profile;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        ProfilePoint p;
        ss >> p.r_cm >> p.rho >> p.temp >> p.ye >> p.pres >> p.entr;

        if (ss) profile.push_back(p);
    }

    return profile;
}

std::vector<SigmaPoint> readSigmaTable(const std::string& filename)
{
    std::ifstream in(filename);
    std::vector<SigmaPoint> table;
    std::string line;

    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        SigmaPoint p;
        ss >> p.Enu_GeV >> p.sigma_GeVminus2 >> p.sigma_cm2;

        if (ss && p.sigma_cm2 > 0.0 && std::isfinite(p.sigma_cm2)) {
            table.push_back(p);
        }
    }

    std::sort(
        table.begin(),
        table.end(),
        [](const SigmaPoint& a, const SigmaPoint& b) {
            return a.Enu_GeV < b.Enu_GeV;
        }
    );

    return table;
}

double interpolateSigmaLogLog(
    const std::vector<SigmaPoint>& table,
    double Enu_GeV
)
{
    if (table.empty()) {
        throw std::runtime_error("Tabela de sigma vazia.");
    }

    if (Enu_GeV <= table.front().Enu_GeV) {
        return table.front().sigma_cm2;
    }

    if (Enu_GeV >= table.back().Enu_GeV) {
        return table.back().sigma_cm2;
    }

    for (size_t i = 0; i + 1 < table.size(); ++i) {
        double E1 = table[i].Enu_GeV;
        double E2 = table[i+1].Enu_GeV;

        if (Enu_GeV >= E1 && Enu_GeV <= E2) {
            double s1 = table[i].sigma_cm2;
            double s2 = table[i+1].sigma_cm2;

            double w =
                (std::log(Enu_GeV) - std::log(E1))
              / (std::log(E2) - std::log(E1));

            return std::exp(std::log(s1) + w*(std::log(s2) - std::log(s1)));
        }
    }

    return table.back().sigma_cm2;
}

std::vector<double> computeTauProfile(
    const std::vector<ProfilePoint>& profile,
    double sigma_cm2
)
{
    int N = static_cast<int>(profile.size());
    std::vector<double> tau(N, 0.0);

    tau[N-1] = 0.0;

    for (int i = N - 2; i >= 0; --i) {
        double r1 = profile[i].r_cm;
        double r2 = profile[i+1].r_cm;

        double rho1 = profile[i].rho;
        double rho2 = profile[i+1].rho;

        double kappa1 = rho1 * NA * sigma_cm2;
        double kappa2 = rho2 * NA * sigma_cm2;

        double dr = r2 - r1;

        tau[i] = tau[i+1] + 0.5*(kappa1 + kappa2)*dr;
    }

    return tau;
}

double findNeutrinosphereRadius(
    const std::vector<ProfilePoint>& profile,
    const std::vector<double>& tau
)
{
    const double tau_sphere = 2.0/3.0;
    int N = static_cast<int>(profile.size());

    for (int i = 0; i < N - 1; ++i) {
        double t1 = tau[i];
        double t2 = tau[i+1];

        if (t1 >= tau_sphere && t2 <= tau_sphere) {
            double r1 = profile[i].r_cm;
            double r2 = profile[i+1].r_cm;

            return r1 + (tau_sphere - t1)*(r2 - r1)/(t2 - t1);
        }
    }

    return -1.0;
}

int main(int argc, char* argv[])
{
    std::string model = "GBW";
    std::string profile_file = "data/collapsar_profile.dat";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--model" && i + 1 < argc) {
            model = argv[++i];
        }
        else if (arg == "--profile" && i + 1 < argc) {
            profile_file = argv[++i];
        }
    }

    if (model != "GBW" && model != "IIM") {
        std::cerr << "Modelo inválido. Use --model GBW ou --model IIM.\n";
        return 1;
    }

    std::string sigma_file =
        "data/sigma_nuN_CC_" + model + ".dat";

    std::string summary_file =
        "data/neutrinosphere_scan_" + model + ".dat";

    std::string profiles_file =
        "data/tau_profiles_scan_" + model + ".dat";

    auto profile = readProfile(profile_file);
    auto sigma_table = readSigmaTable(sigma_file);

    if (profile.empty()) {
        std::cerr << "Erro lendo perfil: " << profile_file << "\n";
        return 1;
    }

    if (sigma_table.empty()) {
        std::cerr << "Erro lendo tabela de sigma: " << sigma_file << "\n";
        return 1;
    }

    std::vector<double> energies = {
        1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11, 1e12
    };

    std::ofstream summary(summary_file);
    summary << "# model " << model << "\n";
    summary << "# profile_file " << profile_file << "\n";
    summary << "# sigma_file " << sigma_file << "\n";
    summary << "# Enu_GeV sigma_cm2 r_nu_cm r_nu_km tau_center\n";

    std::ofstream profiles(profiles_file);
    profiles << "# model " << model << "\n";
    profiles << "# columns: Enu_GeV r_cm r_km tau\n";

    for (double Enu : energies) {
        double sigma_cm2 = interpolateSigmaLogLog(sigma_table, Enu);
        auto tau = computeTauProfile(profile, sigma_cm2);

        double r_nu_cm = findNeutrinosphereRadius(profile, tau);
        double r_nu_km = r_nu_cm > 0.0 ? r_nu_cm/1.0e5 : -1.0;

        summary << Enu << " "
                << sigma_cm2 << " "
                << r_nu_cm << " "
                << r_nu_km << " "
                << tau[0] << "\n";

        for (size_t i = 0; i < profile.size(); ++i) {
            profiles << Enu << " "
                     << profile[i].r_cm << " "
                     << profile[i].r_cm/1.0e5 << " "
                     << tau[i] << "\n";
        }

        profiles << "\n";

        std::cout << "model = " << model
                  << "   E_nu = " << Enu
                  << " GeV, sigma = " << sigma_cm2
                  << " cm^2, r_nu = " << r_nu_km
                  << " km, tau(0) = " << tau[0]
                  << "\n";
    }

    std::cout << "Arquivos gerados:\n";
    std::cout << summary_file << "\n";
    std::cout << profiles_file << "\n";

    return 0;
}