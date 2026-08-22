// =====================================================================
// redshift_scan -- mede
//
//     R(b, E_inf) = [ tau_full - tau_frozen ] / tau_frozen
//
// para duas tabelas de secao de choque total, com a MESMA geometria e a
// MESMA densidade. tau_frozen usa sigma(E_inf); tau_full usa
// sigma(E_loc) = sigma(E_inf/sqrt(f(r))).
//
// POR QUE R, E NAO P NEM tau
//
// O dipolo suprime sigma GLOBALMENTE (saturacao) e ACHATA alpha. Comparar
// tau, ou P = exp(-tau), entre as duas tabelas mistura os dois efeitos, e
// o primeiro domina. Mas
//
//     R = < f^(-alpha/2) > - 1,     peso  n(r) dl
//
// e uma media ponderada de uma razao: um fator global k em sigma
// multiplica numerador e denominador e some. R isola a FORMA.
// Verificado em T39 a 1e-12.
//
// TETO: f(r) >= f(r_t) >= 1/3 ao longo de qualquer raio nao capturado,
// logo R < f(r_t)^(-alpha/2) - 1 <= 3^(alpha/2) - 1, com igualdade so no
// limite b -> b_crit+, onde o raio enrola na esfera de fotons e o peso
// se concentra todo em r_t. Verificado em T37.
// =====================================================================
#include "phasis/cross_section.hpp"
#include "phasis/density.hpp"
#include "phasis/metric.hpp"
#include "phasis/ray.hpp"
#include "phasis/trace.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace phasis;

namespace {

struct Ponto {
    double delta;      // b/b_crit - 1
    double b;
    double r_turn;
    double E;
    double tau_full[2];
    double tau_froz[2];
    double R[2];
    bool   ok[2];
};

} // namespace

int main(int argc, char** argv)
{
    std::string p_dip = "dipole/data/sigma_nuN_CC_GBW.dat";
    std::string p_col = "data/sigma_nuN_CC_collinear.dat";
    std::string out   = "data/redshift_scan.csv";

    double r_s   = 1.0e6;      // cm
    double rho0  = 3.0;        // g/cm^3 em r0
    double r0f   = 1.0e8/1.0e6;// r0 / r_s
    double p_hal = 2.0;
    double r_in_f  = 1.2;      // r_in / r_s  -- ABAIXO da esfera de fotons
    double r_out_f = 1.0e4;    // r_out / r_s
    int    nE = 9;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if      (a == "--dipole")   p_dip   = argv[++i];
        else if (a == "--colinear") p_col   = argv[++i];
        else if (a == "--out")      out     = argv[++i];
        else if (a == "--rs")       r_s     = std::stod(argv[++i]);
        else if (a == "--rho0")     rho0    = std::stod(argv[++i]);
        else if (a == "--r0")       r0f     = std::stod(argv[++i]);
        else if (a == "--p")        p_hal   = std::stod(argv[++i]);
        else if (a == "--rin")      r_in_f  = std::stod(argv[++i]);
        else if (a == "--rout")     r_out_f = std::stod(argv[++i]);
        else if (a == "--nE")       nE      = std::stoi(argv[++i]);
        else { std::cerr << "argumento desconhecido: " << a << "\n"; return 2; }
    }

    try {
        const TableCrossSection dip(p_dip);
        const TableCrossSection col(p_col);
        assert_comparable(dip, col);

        const Schwarzschild schw(r_s);
        const PowerLawHalo halo(rho0, r0f*r_s, p_hal, r_in_f*r_s, r_out_f*r_s);

        const double b_crit = std::sqrt(27.0)/2.0*r_s;

        // A energia LOCAL chega a E/sqrt(f(r_t)) <= sqrt(3) E. A tabela
        // nao extrapola -- por decisao de projeto -- entao o topo da
        // varredura tem de recuar por sqrt(3), senao trace_ray lanca.
        const double E0 = 1.001*std::max(dip.E_min(), col.E_min());
        const double E1 = 0.999*std::min(dip.E_max(), col.E_max())/std::sqrt(3.0);
        if (!(E1 > E0)) { std::cerr << "faixa comum vazia apos recuo sqrt(3)\n"; return 1; }

        // b/b_crit - 1 em escala log: e ai que o efeito cresce.
        std::vector<double> deltas;
        for (int k = 0; k <= 10; ++k) deltas.push_back(std::pow(10.0, -12.0 + 1.2*k));

        IntegratorOpts opt;

        const CrossSection* xs[2] = { &dip, &col };
        std::vector<Ponto> pts;

        for (double d : deltas) {
            const double b = b_crit*(1.0 + d);
            for (int i = 0; i < nE; ++i) {
                const double t = (nE == 1) ? 0.0 : double(i)/(nE - 1);
                Ponto q;
                q.delta = d;
                q.b     = b;
                q.E     = E0*std::pow(E1/E0, t);
                q.r_turn = schw.r_turning(b, b);
                for (int m = 0; m < 2; ++m) {
                    Ray ray; ray.E_inf_GeV = q.E; ray.b_cm = b;
                    const Result a = trace_ray(ray, schw, halo, *xs[m], opt, false);
                    const Result c = trace_ray(ray, schw, halo, *xs[m], opt, true);
                    q.tau_full[m] = a.tau;
                    q.tau_froz[m] = c.tau;
                    q.R[m] = (c.tau > 0.0) ? (a.tau - c.tau)/c.tau : 0.0;
                    q.ok[m] = a.tolerance_met && c.tolerance_met && !a.captured;
                }
                pts.push_back(q);
            }
        }

        std::ofstream f(out);
        if (!f) { std::cerr << "nao consegui escrever " << out << "\n"; return 1; }
        f.precision(12);
        f << "# redshift_scan -- R(b,E) = (tau_full - tau_frozen)/tau_frozen\n";
        f << "# mesma geometria e mesma densidade nos dois; so o argumento\n";
        f << "# de sigma muda (E_loc contra E_inf).\n";
        f << "# metrica = Schwarzschild  r_s_cm = " << r_s
          << "  b_crit_cm = " << b_crit << "\n";
        f << "# perfil = PowerLawHalo  rho0 = " << rho0
          << "  r0_cm = " << r0f*r_s << "  p = " << p_hal
          << "  r_in_cm = " << r_in_f*r_s << "  r_out_cm = " << r_out_f*r_s << "\n";
        f << "# nota: r_in = " << r_in_f << " r_s fica ABAIXO da esfera de fotons\n";
        f << "#       (1.5 r_s), sem o que os raios quase-criticos enrolariam\n";
        f << "#       no vazio e o teto 3^(alpha/2) seria inatingivel.\n";
        write_metadata_header(f, "dipole",   dip);
        write_metadata_header(f, "colinear", col);
        f << "b_over_bcrit_minus_1,b_cm,r_turn_over_rs,E_inf_GeV,"
             "tau_dipole_full,tau_dipole_frozen,R_dipole,"
             "tau_colinear_full,tau_colinear_frozen,R_colinear,ok\n";
        for (const auto& q : pts) {
            f << q.delta << "," << q.b << "," << q.r_turn/r_s << "," << q.E << ","
              << q.tau_full[0] << "," << q.tau_froz[0] << "," << q.R[0] << ","
              << q.tau_full[1] << "," << q.tau_froz[1] << "," << q.R[1] << ","
              << ((q.ok[0] && q.ok[1]) ? 1 : 0) << "\n";
        }

        std::printf("%-11s %-11s %-12s %-11s %-11s %-9s %s\n",
                    "b/bc - 1", "r_t/r_s", "E_inf", "R_dipolo", "R_colinear",
                    "a_dip", "a_col");
        for (const auto& q : pts) {
            if (q.E < 0.9e7 || q.E > 1.1e7) continue;
            std::printf("%-11.1e %-11.5f %-12.4g %-11.5f %-11.5f %-9.4f %.4f\n",
                        q.delta, q.r_turn/r_s, q.E, q.R[0], q.R[1],
                        dip.log_slope(q.E), col.log_slope(q.E));
        }
        std::cout << "\nArquivo gerado: " << out << "\n";
    } catch (const std::exception& e) {
        std::cerr << "ERRO: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
