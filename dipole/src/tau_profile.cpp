#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>

constexpr double NA = 6.02214076e23;

struct ProfilePoint {
    double r_cm;
    double rho;
    double temp;
    double ye;
    double pres;
    double entr;
};

struct SigmaPoint {
    double Enu_GeV;
    double sigma_GeVminus2;
    double sigma_cm2;
};

struct TauPoint {
    double r_cm;
    double rho;
    double tau;
    double survival;
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

        if (ss) table.push_back(p);
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
        const double E1 = table[i].Enu_GeV;
        const double E2 = table[i+1].Enu_GeV;

        if (Enu_GeV >= E1 && Enu_GeV <= E2) {
            const double s1 = table[i].sigma_cm2;
            const double s2 = table[i+1].sigma_cm2;

            const double logE  = std::log(Enu_GeV);
            const double logE1 = std::log(E1);
            const double logE2 = std::log(E2);

            const double logs1 = std::log(s1);
            const double logs2 = std::log(s2);

            const double w = (logE - logE1)/(logE2 - logE1);

            return std::exp(logs1 + w*(logs2 - logs1));
        }
    }

    return table.back().sigma_cm2;
}

int main(int argc, char* argv[])
{
    std::string profile_file = "data/flash_like_profile.dat";
    std::string sigma_file   = "data/sigma_nuN_CC_GBW.dat";
    std::string output_file  = "data/tau_profile.dat";

    double Enu_GeV = 1.0e9;

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
        else if (arg == "--Enu" && i + 1 < argc) {
            Enu_GeV = std::stod(argv[++i]);
        }
    }

    auto profile = readProfile(profile_file);
    auto sigma_table = readSigmaTable(sigma_file);

    if (profile.empty()) {
        std::cerr << "Erro lendo perfil: " << profile_file << "\n";
        return 1;
    }

    if (sigma_table.empty()) {
        std::cerr << "Erro lendo tabela sigma: " << sigma_file << "\n";
        return 1;
    }

    const double sigma_cm2 =
        interpolateSigmaLogLog(sigma_table, Enu_GeV);

    const int N = static_cast<int>(profile.size());
    std::vector<TauPoint> tau_profile(N);

    tau_profile[N-1].r_cm = profile[N-1].r_cm;
    tau_profile[N-1].rho = profile[N-1].rho;
    tau_profile[N-1].tau = 0.0;
    tau_profile[N-1].survival = 1.0;

    for (int i = N - 2; i >= 0; --i) {
        const double r1 = profile[i].r_cm;
        const double r2 = profile[i+1].r_cm;

        const double rho1 = profile[i].rho;
        const double rho2 = profile[i+1].rho;

        const double kappa1 = rho1 * NA * sigma_cm2;
        const double kappa2 = rho2 * NA * sigma_cm2;

        const double dr = r2 - r1;

        const double dtau =
            0.5*(kappa1 + kappa2)*dr;

        tau_profile[i].r_cm = r1;
        tau_profile[i].rho = rho1;
        tau_profile[i].tau = tau_profile[i+1].tau + dtau;

        if (tau_profile[i].tau < 700.0) {
            tau_profile[i].survival =
                std::exp(-tau_profile[i].tau);
        } else {
            tau_profile[i].survival = 0.0;
        }
    }

    std::ofstream out(output_file);

    out << "# profile_file " << profile_file << "\n";
    out << "# sigma_table " << sigma_file << "\n";
    out << "# Enu_GeV " << Enu_GeV << "\n";
    out << "# sigma_cm2_interpolated " << sigma_cm2 << "\n";
    out << "# r_cm rho_gcm3 tau survival\n";

    for (const auto& p : tau_profile) {
        out << p.r_cm << " "
            << p.rho << " "
            << p.tau << " "
            << p.survival << "\n";
    }

    double r_nu = -1.0;
    const double tau_sphere = 2.0/3.0;

    for (int i = 0; i < N - 1; ++i) {
        const double t1 = tau_profile[i].tau;
        const double t2 = tau_profile[i+1].tau;

        if (t1 >= tau_sphere && t2 <= tau_sphere) {
            const double r1 = tau_profile[i].r_cm;
            const double r2 = tau_profile[i+1].r_cm;

            r_nu =
                r1 + (tau_sphere - t1)*(r2 - r1)/(t2 - t1);

            break;
        }
    }

    std::cout << "Perfil usado: " << profile_file << "\n";
    std::cout << "Tabela sigma usada: " << sigma_file << "\n";
    std::cout << "E_nu = " << Enu_GeV << " GeV\n";
    std::cout << "sigma_nuN interpolado = " << sigma_cm2 << " cm^2\n";
    std::cout << "Arquivo gerado: " << output_file << "\n";

    if (r_nu > 0.0) {
        std::cout << "Neutrinosfera aproximada:\n";
        std::cout << "r_nu = " << r_nu/1.0e5 << " km\n";
    } else {
        std::cout << "Neutrinosfera nao encontrada.\n";
    }

    return 0;
}