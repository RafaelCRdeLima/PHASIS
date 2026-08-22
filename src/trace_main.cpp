// =====================================================================
// Executavel `trace`.
//
//   trace config.txt [saida.csv]
//
// Le um arquivo chave=valor e escreve CSV com
//
//   E_inf_GeV,b_cm,r_min_cm,tau,P_surv,column_g_cm2,path_cm,crosses_matter
//
// Chaves reconhecidas estao listadas em `imprime_ajuda()`.
// =====================================================================

#include "phasis/cross_section.hpp"
#include "phasis/density.hpp"
#include "phasis/metric.hpp"
#include "phasis/trace.hpp"
#include "phasis/units.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace phasis;

namespace {

using Config = std::map<std::string, std::string>;

Config le_config(const std::string& caminho)
{
    std::ifstream in(caminho);
    if (!in) throw std::runtime_error("nao consegui abrir " + caminho);

    Config c;
    std::string linha;
    while (std::getline(in, linha)) {
        const auto h = linha.find('#');
        if (h != std::string::npos) linha = linha.substr(0, h);

        const auto eq = linha.find('=');
        if (eq == std::string::npos) continue;

        auto trim = [](std::string s) {
            const auto a = s.find_first_not_of(" \t\r\n");
            if (a == std::string::npos) return std::string{};
            const auto b = s.find_last_not_of(" \t\r\n");
            return s.substr(a, b - a + 1);
        };

        const std::string k = trim(linha.substr(0, eq));
        const std::string v = trim(linha.substr(eq + 1));
        if (!k.empty()) c[k] = v;
    }
    return c;
}

std::string get(const Config& c, const std::string& k, const std::string& def)
{
    const auto it = c.find(k);
    return (it == c.end()) ? def : it->second;
}

double getd(const Config& c, const std::string& k, double def)
{
    const auto it = c.find(k);
    if (it == c.end()) return def;
    return std::stod(it->second);
}

void imprime_ajuda()
{
    std::cout <<
R"(uso: trace CONFIG [SAIDA.csv]

Chaves do arquivo de configuracao (chave = valor, '#' comenta):

  metric        = minkowski | schwarzschild
    M_solar     = 10.0                      (schwarzschild) OU
    r_s_cm      = 2.95325e6                 (schwarzschild)

  density       = uniform_ball | powerlaw_halo
    rho0_g_cm3  = 5.51
    R_cm        = 6.371e8                   (uniform_ball)
    r0_cm       = 1e8                       (powerlaw_halo)
    p           = 2.0                       (powerlaw_halo)
    r_in_cm     = 1e7                       (powerlaw_halo)
    r_out_cm    = 1e9                       (powerlaw_halo)

  xsec          = powerlaw | table
    sigma0_cm2  = 1e-33                     (powerlaw)
    E0_GeV      = 1e3                       (powerlaw)
    alpha       = 0.363                     (powerlaw)
    table_path  = dipole/data/sigma_nuN_CC_GBW.dat      (table)

  xsec_current  = cc | total
    As tabelas do dipole sao CC PURAS. Para 'total' e obrigatorio dar
    UM dos dois abaixo -- nao ha razao NC/CC assumida por default:
    table_path_nc  = ...                    tabela NC separada, OU
    nc_to_cc_ratio = 0.42                   razao declarada explicitamente

  E_inf_GeV     = 1e9                       energia unica, OU:
  E_min_GeV     = 1e4
  E_max_GeV     = 1e12
  n_E           = 25                        (log-espacado)

  b_cm          = 0.0                       impacto unico, OU:
  b_min_cm      = 0.0
  b_max_cm      = 6.371e8
  n_b           = 20                        (linear)

  rel_tol       = 1e-8
  max_depth     = 50
  max_evals     = 2e6                       orcamento global da quadratura
  deflection    = 1                          0 desliga o calculo do angulo
)";
}

std::vector<double> grade(const Config& c,
                          const std::string& unico,
                          const std::string& kmin,
                          const std::string& kmax,
                          const std::string& kn,
                          bool log_espacado,
                          double def)
{
    if (c.count(unico)) return { getd(c, unico, def) };

    const int n = static_cast<int>(getd(c, kn, 1));
    const double a = getd(c, kmin, def);
    const double b = getd(c, kmax, def);

    std::vector<double> v;
    if (n <= 1) { v.push_back(a); return v; }

    v.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double t = static_cast<double>(i)/(n - 1);
        v.push_back(log_espacado ? a*std::pow(b/a, t) : a + t*(b - a));
    }
    return v;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) { imprime_ajuda(); return 2; }

    try {
        const Config cfg = le_config(argv[1]);

        // ---- metrica --------------------------------------------------
        std::unique_ptr<Metric> metric;
        const std::string mname = get(cfg, "metric", "minkowski");
        if (mname == "minkowski") {
            metric = std::make_unique<Minkowski>();
        } else if (mname == "schwarzschild") {
            if (cfg.count("M_solar")) {
                metric = std::make_unique<Schwarzschild>(
                    Schwarzschild::from_solar_masses(getd(cfg, "M_solar", 1.0)));
            } else {
                metric = std::make_unique<Schwarzschild>(getd(cfg, "r_s_cm", 1.0e5));
            }
        } else {
            throw std::runtime_error("metric desconhecida: " + mname);
        }

        // ---- densidade ------------------------------------------------
        std::unique_ptr<DensityProfile> profile;
        const std::string dname = get(cfg, "density", "uniform_ball");
        if (dname == "uniform_ball") {
            profile = std::make_unique<UniformBall>(
                getd(cfg, "rho0_g_cm3", units::earth_mean_density),
                getd(cfg, "R_cm", units::earth_radius_cm));
        } else if (dname == "powerlaw_halo") {
            profile = std::make_unique<PowerLawHalo>(
                getd(cfg, "rho0_g_cm3", 1.0),
                getd(cfg, "r0_cm", 1.0e8),
                getd(cfg, "p", 2.0),
                getd(cfg, "r_in_cm", 1.0e7),
                getd(cfg, "r_out_cm", 1.0e9));
        } else {
            throw std::runtime_error("density desconhecida: " + dname);
        }

        // ---- secao de choque ------------------------------------------
        std::shared_ptr<const CrossSection> xsec;
        const std::string xname = get(cfg, "xsec", "powerlaw");
        if (xname == "powerlaw") {
            xsec = std::make_shared<PowerLawCrossSection>(
                getd(cfg, "sigma0_cm2", 1.0e-33),
                getd(cfg, "E0_GeV", 1.0e3),
                getd(cfg, "alpha", 0.363));
        } else if (xname == "table") {
            const std::string p = get(cfg, "table_path", "");
            if (p.empty()) throw std::runtime_error("xsec=table exige table_path");
            xsec = std::make_shared<TableCrossSection>(p);
        } else {
            throw std::runtime_error("xsec desconhecida: " + xname);
        }

        // Corrente carregada apenas, ou total?
        //
        // As tabelas de PHASIS/dipole sao CC PURAS. Na Fase 1/2, sem
        // regeneracao, as duas convencoes sao defensaveis por motivos
        // OPOSTOS: sigma_CC porque so a corrente carregada remove o
        // neutrino de vez; sigma_tot porque, sem regeneracao, a corrente
        // neutra tambem tira o neutrino do bin de energia.
        //
        // Por isso a escolha e explicita e vai para o cabecalho do CSV.
        // Nao ha default silencioso.
        const std::string corrente = get(cfg, "xsec_current", "cc");
        std::string nota_corrente = "cc (tabela como esta)";

        if (corrente == "total") {
            const std::string p_nc = get(cfg, "table_path_nc", "");
            if (!p_nc.empty()) {
                auto soma = std::make_shared<SumCrossSection>();
                soma->add(xsec);
                soma->add(std::make_shared<TableCrossSection>(p_nc));
                xsec = soma;
                nota_corrente = "total = CC + NC (tabela NC: " + p_nc + ")";
            } else if (cfg.count("nc_to_cc_ratio")) {
                const double k = 1.0 + getd(cfg, "nc_to_cc_ratio", 0.0);
                xsec = std::make_shared<ScaledCrossSection>(xsec, k);
                std::ostringstream m;
                m << "total = (1 + " << (k - 1.0) << ") * CC  [razao declarada]";
                nota_corrente = m.str();
            } else {
                throw std::runtime_error(
                    "xsec_current=total exige table_path_nc OU nc_to_cc_ratio "
                    "explicito. As tabelas do dipole sao CC puras, e nao ha "
                    "razao NC/CC assumida por default.");
            }
        } else if (corrente != "cc") {
            throw std::runtime_error("xsec_current deve ser cc ou total");
        }

        // ---- grades ---------------------------------------------------
        const auto Es = grade(cfg, "E_inf_GeV", "E_min_GeV", "E_max_GeV", "n_E", true, 1.0e9);
        const auto bs = grade(cfg, "b_cm", "b_min_cm", "b_max_cm", "n_b", false, 0.0);

        IntegratorOpts opts;
        opts.rel_tol   = getd(cfg, "rel_tol", 1.0e-8);
        opts.max_depth = static_cast<int>(getd(cfg, "max_depth", 50));
        opts.max_evals = static_cast<long>(getd(cfg, "max_evals", 2.0e6));
        opts.want_deflection = (getd(cfg, "deflection", 1.0) != 0.0);

        // ---- saida ----------------------------------------------------
        std::ostream* out = &std::cout;
        std::ofstream arq;
        if (argc >= 3) {
            arq.open(argv[2]);
            if (!arq) throw std::runtime_error(std::string("nao consegui escrever ") + argv[2]);
            out = &arq;
        }

        *out << "# PHASIS trace\n";
        *out << "# metric=" << metric->name()
             << " density=" << profile->name()
             << " xsec=" << xsec->name() << "\n";
        *out << "# corrente=" << nota_corrente << "\n";
        if (const auto* sw = dynamic_cast<const Schwarzschild*>(metric.get())) {
            *out << "# r_s=" << sw->r_s() << " cm  b_crit=" << sw->b_crit()
                 << " cm  r_photon=" << sw->r_photon() << " cm\n";
        }
        *out << "# unidades: cm, g/cm^3, cm^2, GeV\n";
        *out << "E_inf_GeV,b_cm,r_min_cm,tau,P_surv,column_g_cm2,path_cm,"
                "crosses_matter,captured,near_critical,deflection_rad,"
                "E_loc_max_GeV,winding_turns,tolerance_met\n";
        out->precision(12);

        for (double E : Es) {
            for (double b : bs) {
                Ray ray; ray.E_inf_GeV = E; ray.b_cm = b;
                const Result r = trace_ray(ray, *metric, *profile, *xsec, opts);
                *out << E << "," << b << "," << r.r_min_cm << ","
                     << r.tau << "," << r.P_surv << ","
                     << r.column_density << "," << r.path_length_cm << ","
                     << (r.crosses_matter ? 1 : 0) << ","
                     << (r.captured ? 1 : 0) << ","
                     << (r.near_critical ? 1 : 0) << ","
                     << r.deflection_rad << "," << r.E_loc_max_GeV << ","
                     << r.winding_turns << ","
                     << (r.tolerance_met ? 1 : 0) << "\n";
            }
        }

        if (argc >= 3) std::cerr << "escrito: " << argv[2] << "\n";
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "erro: " << e.what() << "\n";
        return 1;
    }
}
