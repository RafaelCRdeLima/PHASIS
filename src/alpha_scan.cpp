// =====================================================================
// alpha_scan -- tabula alpha_eff(E) = d ln sigma / d ln E e o teto
// 3^(alpha/2) para duas tabelas de secao de choque TOTAL.
//
// A grandeza medida e propriedade so de sigma_tot: nao depende de NC,
// nem de dsigma/dy, nem do acoplamento NC do dipole. sigma_CC ja esta
// validada, entao esta medida esta desbloqueada hoje.
//
// O teto vem de E_loc/E_inf = 1/sqrt(f(r_t)) <= sqrt(3) em
// Schwarzschild: todo raio nao capturado tem r_t >= 1.5 r_s, onde
// f = 1/3. Com sigma ~ E^alpha, a razao de secoes de choque entre o
// ponto de retorno e o infinito e no maximo 3^(alpha/2).
//
// A ASSINATURA nao e o valor do teto, e a DERIVADA dele em E. Quanto
// mais fundo na saturacao, mais achatado alpha, e mais baixo o teto.
// Um deslocamento de normalizacao qualquer reajuste de PDF reproduz;
// uma inclinacao que cai com a energia, nao.
// =====================================================================
#include "phasis/cross_section.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace phasis;

int main(int argc, char** argv)
{
    std::string p_dip = "dipole/data/sigma_nuN_CC_GBW.dat";
    std::string p_col = "data/sigma_nuN_CC_collinear.dat";
    std::string out   = "data/alpha_scan.csv";
    int nE = 120;

    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if      (a == "--dipole")   p_dip = argv[++i];
        else if (a == "--colinear") p_col = argv[++i];
        else if (a == "--out")      out   = argv[++i];
        else if (a == "--nE")       nE    = std::stoi(argv[++i]);
        else { std::cerr << "argumento desconhecido: " << a << "\n"; return 2; }
    }

    try {
        const TableCrossSection dip(p_dip);
        const TableCrossSection col(p_col);

        // O UNICO ponto onde uma comparacao errada passaria despercebida.
        assert_comparable(dip, col);

        // Faixa COMUM as duas. Fora dela nenhuma das duas extrapola --
        // e nao se compara o que uma so delas sabe.
        const double E0 = std::max(dip.E_min(), col.E_min());
        const double E1 = std::min(dip.E_max(), col.E_max());
        if (!(E1 > E0)) {
            std::cerr << "as duas tabelas nao tem faixa em comum\n";
            return 1;
        }

        std::ofstream f(out);
        if (!f) { std::cerr << "nao consegui escrever " << out << "\n"; return 1; }
        f.precision(10);

        f << "# alpha_scan -- alpha_eff(E) = d ln sigma / d ln E e teto 3^(alpha/2)\n";
        f << "# estimador: diferenca centrada em ln E, passo adaptativo por\n";
        f << "#            reducao pela metade (h inicial 0.1, piso 3e-3).\n";
        f << "#            A rota analitica esta fechada: a interpolacao da\n";
        f << "#            tabela e linear em log-log, portanto C^0, e a\n";
        f << "#            derivada dela e uma funcao escada.\n";
        f << "# faixa_comum_GeV = " << E0 << " " << E1 << "\n";
        write_metadata_header(f, "dipole",   dip);
        write_metadata_header(f, "colinear", col);

        int nao_conv = 0;
        std::vector<std::string> linhas;
        for (int i = 0; i < nE; ++i) {
            const double t = (nE == 1) ? 0.0 : double(i)/(nE - 1);
            const double E = E0*std::pow(E1/E0, t);

            const LogSlope ad = dip.log_slope_detail(E);
            const LogSlope ac = col.log_slope_detail(E);
            if (!ad.converged || !ac.converged) ++nao_conv;

            char buf[256];
            std::snprintf(buf, sizeof buf, "%.10g,%.10g,%.10g,%.10g,%.10g",
                          E, ad.alpha, ac.alpha,
                          std::pow(3.0, 0.5*ad.alpha),
                          std::pow(3.0, 0.5*ac.alpha));
            linhas.emplace_back(buf);
        }

        f << "# pontos_sem_convergencia_do_passo = " << nao_conv
          << " de " << nE << "\n";
        f << "E_GeV,alpha_dipole,alpha_colinear,teto_dipole,teto_colinear\n";
        for (const auto& l : linhas) f << l << "\n";

        std::printf("%-13s %-11s %-13s %-11s %-13s\n",
                    "E_GeV", "a_dipolo", "a_colinear", "teto_dip", "teto_col");
        for (int i = 0; i < nE; i += std::max(1, nE/12)) {
            const double t = (nE == 1) ? 0.0 : double(i)/(nE - 1);
            const double E = E0*std::pow(E1/E0, t);
            const double a = dip.log_slope(E), c = col.log_slope(E);
            std::printf("%-13.4g %-11.5f %-13.5f %-11.5f %-13.5f\n",
                        E, a, c, std::pow(3.0,0.5*a), std::pow(3.0,0.5*c));
        }
        std::cout << "\nArquivo gerado: " << out << "\n";
    } catch (const std::exception& e) {
        std::cerr << "ERRO: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
