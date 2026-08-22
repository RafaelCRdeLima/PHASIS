#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "parameters.hpp"
#include "dipole_models.hpp"
#include "integrals.hpp"
#include "structure_table.hpp"
#include "sigma_nuN_core.hpp"
#include "weak_structure_functions.hpp"

using namespace dipole;

enum class DipoleModel { GBW, IIM };

static DipoleModel parseDipoleModel(const std::string& model_name)
{
    if (model_name == "GBW") return DipoleModel::GBW;
    if (model_name == "IIM") return DipoleModel::IIM;
    throw std::runtime_error("Modelo invalido. Use --model GBW ou --model IIM.");
}

int main(int argc, char* argv[])
{
    bool useF3 = false;
    std::string pdf_set = "CT10nlo";
    std::string beam_name = "nu";
    double logEmin = 3.0;
    double logEmax = 14.0;
    int NE = 20;

    // F4: densidade de nos por decada, em vez de numero fixo.
    double nodesQ = 16.0;
    double nodesX = 16.0;

    // Quadratura interna em (r,z) usada para montar a tabela.
    QuadratureGrid quad;
    quad.Nr = 120;
    quad.Nz = 120;

    TableSpec spec;

    double Q2min = 1.0;
    std::string model_name = "GBW";
    std::string current_name = "CC";
    int nY = 0;                 // > 0: grava tambem a tabela dsigma/dy
    double yMinTable = 0.0;     // 0 = deriva do corte cinematico
    bool verbose = true;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if      (arg == "--model")    model_name = argv[++i];
        else if (arg == "--logEmin")  logEmin = std::stod(argv[++i]);
        else if (arg == "--logEmax")  logEmax = std::stod(argv[++i]);
        else if (arg == "--NE")       NE = std::stoi(argv[++i]);
        else if (arg == "--nodesQ")   nodesQ = std::stod(argv[++i]);
        else if (arg == "--nodesX")   nodesX = std::stod(argv[++i]);
        else if (arg == "--Nr")       quad.Nr = std::stoi(argv[++i]);
        else if (arg == "--Nz")       quad.Nz = std::stoi(argv[++i]);
        else if (arg == "--Ntabx")    spec.Nx = std::stoi(argv[++i]);
        else if (arg == "--NtabQ")    spec.NQ = std::stoi(argv[++i]);
        else if (arg == "--Q2min")    Q2min = std::stod(argv[++i]);
        else if (arg == "--use-F3")   useF3 = std::stoi(argv[++i]) != 0;
        else if (arg == "--pdf-set")  pdf_set = argv[++i];
        else if (arg == "--beam")     beam_name = argv[++i];
        else if (arg == "--current")  current_name = argv[++i];
        else if (arg == "--nY")       nY = std::stoi(argv[++i]);
        else if (arg == "--yMin")     yMinTable = std::stod(argv[++i]);
        else if (arg == "--quiet")    verbose = false;
        else {
            std::cerr << "Argumento desconhecido: " << arg << "\n";
            return 2;
        }
    }

    if (NE < 1) {
        std::cerr << "--NE precisa ser >= 1.\n";
        return 2;
    }

    DipoleModel model = parseDipoleModel(model_name);

    CurrentType current;
    if      (current_name == "CC") current = CurrentType::CC;
    else if (current_name == "NC") current = CurrentType::NC;
    else { std::cerr << "Corrente invalida. Use --current CC ou --current NC.\n"; return 2; }

    GBWParameters gbw;
    IIMParameters iim;

    // As massas efetivas fazem parte do ajuste de cada modelo.
    QuarkMasses masses;
    if (model == DipoleModel::GBW) {
        masses.u = masses.d = masses.s = gbw.m_light;
        masses.c = gbw.m_charm;
    } else {
        masses.u = masses.d = masses.s = iim.m_light;
        masses.c = iim.m_charm;
    }

    weak::BeamType beam = weak::parseBeamType(beam_name);

    std::unique_ptr<weak::WeakStructureFunctions> weakSF;
    if (useF3) {
        weakSF = std::make_unique<weak::WeakStructureFunctions>(pdf_set, 0);
    }

    // A tabela precisa cobrir todo o (x, Q^2) que o laco de energia
    // vai pedir: x_min = Q2min/s_max e Q2_max = 0.999*s_max.
    const double s_max = 2.0*MN*std::pow(10.0, logEmax);
    spec.Q2Min = Q2min;
    spec.Q2Max = 1.05*0.999*s_max;
    spec.xMin  = 0.5*Q2min/s_max;
    spec.xMax  = 1.0;

    // CC -- Kutak-Kwiecinski: os dipolos favorecidos por Cabibbo sao
    // ud~(du~) e cs~(sc~). Os dois canais somados dao o Sigma(carga^2)
    // efetivo = 4 das eqs. (12) e (13).
    //
    // NC -- o Z acopla ao MESMO sabor dos dois lados, entao sao quatro
    // dipolos qq~ em vez de dois. A soma efetiva vira
    // Soma_q (g_V^2+g_A^2) = 1.3127, que e o C de KK eq. (15). A razao
    // 1.3127/4 = 0.328, vezes (M_Z/M_W)^2, da 0.423 -- a razao
    // sigma_NC/sigma_CC conhecida. Ver o bloco D10 em parameters.hpp.
    //
    // Mesmo conjunto de sabores nos dois casos (u, d, s, c): o ajuste
    // GBW e de quatro sabores, e incluir b seria sair do ajuste.
    const std::vector<StructureTable::Channel> channels_cc = {
        { QuarkFlavor::u, QuarkFlavor::d },
        { QuarkFlavor::c, QuarkFlavor::s }
    };
    const std::vector<StructureTable::Channel> channels_nc = {
        { QuarkFlavor::u, QuarkFlavor::u },
        { QuarkFlavor::d, QuarkFlavor::d },
        { QuarkFlavor::s, QuarkFlavor::s },
        { QuarkFlavor::c, QuarkFlavor::c }
    };
    const std::vector<StructureTable::Channel>& channels =
        (current == CurrentType::NC) ? channels_nc : channels_cc;

    if (verbose) {
        std::cerr << "Montando a tabela de F_T, F_L (uma vez, serve todas as energias)\n";
    }

    const StructureTable table(
        model == DipoleModel::GBW ? DipoleModelId::GBW : DipoleModelId::IIM,
        channels, gbw, iim, quad, spec, masses, verbose, current
    );

    if (verbose) {
        std::cerr << "  " << table.describe() << "\n"
                  << "  construida em " << table.buildSeconds() << " s\n\n";
    }

    const std::string output_file =
        "data/sigma_nuN_" + current_name + "_" + model_name + ".dat";

    std::ofstream out(output_file);
    out << std::setprecision(10);

    // F9: proveniencia. A tabela tem de poder reproduzir a si mesma
    // so pelo cabecalho.
    //
    // As NOVE chaves obrigatorias vem primeiro, no formato
    // '# chave = valor'. Elas nao sao decoracao: sao o que
    // phasis::assert_comparable exige antes de deixar comparar esta
    // tabela com outra. Uma tabela de nu contra uma de nubar, ou de
    // alvo isoescalar contra proton, continuam sendo secoes de choque
    // validas e a razao entre elas continua sendo um numero -- so a
    // comparacao e que deixa de significar algo, e nada DENTRO de cada
    // arquivo consegue detectar isso.
    out << "# convention_y = (E_in - E_out)/E_in\n";
    out << "# target       = isoscalar_nucleon\n";
    out << "# projectile   = " << beam_name << "\n";
    out << "# current      = " << current_name << "\n";
    out << "# units_sigma  = cm^2\n";
    out << "# units_E      = GeV\n";
    out << "# M_Z_GeV      = " << MZ << "\n";
    out << "# dipole_model = " << model_name << "\n";
    out << "# generated_by = sigma_nuN (PHASIS/dipole) commit "
        << PHASIS_GIT_HASH << "\n";
    // extras, fora das nove
    out << "# M_W_GeV      = " << MW << "\n";
    out << "# xF3_from     = " << (useF3 ? pdf_set : std::string("nenhum")) << "\n";
    out << "# model " << model_name << "\n";
    out << "# gerado_por sigma_nuN (PHASIS/dipole)\n";
    out << "# formalismo Kutak-Kwiecinski EPJ C29 (2003) 521, eqs (2),(8),(9),(12),(13)\n";
    out << "# beam " << beam_name << "\n";
    out << "# use_F3 " << (useF3 ? 1 : 0);
    if (useF3) out << "  pdf_set " << pdf_set;
    out << "\n";
    out << "# Q2min_GeV2 " << Q2min << "\n";
    out << "# largeXFactor (1-x)^7   [regra de contagem de constituintes, n_s=4]\n";
    if (model == DipoleModel::GBW) {
        out << "# GBW [Golec-Biernat & Wusthoff, PRD 59 (1999) 014017; ajuste de 4 sabores]\n";
        out << "# GBW sigma0_mb " << gbw.sigma0_mb
            << "  lambda " << gbw.lambda
            << "  x0 " << gbw.x0
            << "  Q0sq " << gbw.Q0sq << "\n";
    } else {
        out << "# bCGC [Rezaeian & Schmidt, PRD 88 (2013) 074016, Tabela II, m_c=1.4]\n";
        out << "# bCGC B_CGC " << iim.BCGC
            << "  gamma_s " << iim.gamma_s
            << "  N0 " << iim.N0
            << "  x0 " << iim.x0
            << "  lambda " << iim.lambda
            << "  bMax " << iim.bMax
            << "  Nb " << iim.Nb << "\n";
    }
    out << "# massas_GeV  u " << masses.u << "  d " << masses.d
        << "  s " << masses.s << "  c " << masses.c << "\n";
    out << "# canais " << (current == CurrentType::NC ? "uu, dd, ss, cc"
                                                       : "ud, cs") << "\n";
    out << "# quadratura_rz  Nr " << quad.Nr << "  Nz " << quad.Nz
        << "  rMin " << quad.rMin << "  rMax " << quad.rMax
        << "  zMin " << quad.zMin << "   [ln r; z log nas duas pontas]\n";
    out << "# tabela_F  Nx " << spec.Nx << "  NQ " << spec.NQ
        << "  xMin " << spec.xMin << "  xMax " << spec.xMax
        << "  Q2Min " << spec.Q2Min << "  Q2Max " << spec.Q2Max << "\n";
    out << "# integracao_sigma  nos_por_decada_Q " << nodesQ
        << "  nos_por_decada_x " << nodesX << "\n";
    out << "# Enu_GeV sigma_GeV_minus2 sigma_cm2\n";

    for (int i = 0; i < NE; ++i) {

        const double logE = (NE == 1)
            ? logEmin
            : logEmin + i*(logEmax - logEmin)/(NE - 1);

        const double Enu = std::pow(10.0, logE);

        int nodesQused = 0;
        const double sigma_gev2 = sigmaNuN(
            Enu, current, table, nodesQ, nodesX, useF3, weakSF.get(), beam,
            Q2min, &nodesQused
        );

        const double sigma_cm2 = sigma_gev2*GeVminus2_to_cm2;

        out << Enu << " " << sigma_gev2 << " " << sigma_cm2 << "\n";

        if (verbose) {
            std::cout << "model = " << model_name
                      << "   E_nu = " << Enu
                      << " GeV   sigma_" << current_name << " = " << sigma_cm2
                      << " cm^2   (NlogQ = " << nodesQused << ")\n";
        }
    }

    // =================================================================
    // Tabela diferencial dsigma/dy, no formato que
    // phasis::TableDifferentialCrossSection consome.
    //
    // GRADE EM y, e o que ela NAO cobre
    //
    // dsigma/dy so e nao nula para y >= Q2min/(s x_max): abaixo disso
    // todo o Q^2 acessivel cai sob o corte de validade do modelo. Esse
    // limite DEPENDE DA ENERGIA (cai como 1/E), mas o formato de tabela
    // exige uma grade em y comum a todas as energias.
    //
    // A escolha aqui e cobrir tudo: y_min da grade e o corte na energia
    // MAIS ALTA, e nas energias mais baixas as linhas abaixo do proprio
    // corte saem zeradas. Sao zeros de verdade -- consequencia de
    // Q^2 >= Q2min -- e nao lacunas. A alternativa (comecar a grade no
    // corte da energia mais BAIXA) esconderia, nas energias altas, uma
    // faixa de y pequeno que o modelo cobre e que a cascata usa: y
    // pequeno e perda de energia pequena, que e o regime onde o kernel
    // de regeneracao mais pesa.
    //
    // O cabecalho reporta, por energia, que fracao de sigma a integral
    // em y da grade recupera. Se essa fracao nao for ~1, o problema e a
    // RESOLUCAO da grade, e o numero esta la para ser visto.
    // =================================================================
    if (nY > 1) {
        const double E_hi = std::pow(10.0, logEmax);
        const double y_lo = (yMinTable > 0.0)
            ? yMinTable
            : yMinKinematic(E_hi, Q2min);

        const std::string dy_file =
            "data/dsigma_dy_" + current_name + "_" + model_name + ".dat";

        std::ofstream dy(dy_file);
        if (!dy) {
            std::cerr << "Nao consegui escrever " << dy_file << "\n";
            return 1;
        }
        dy << std::setprecision(10);

        dy << "# convention_y = (E_in - E_out)/E_in\n";
        dy << "# target       = isoscalar_nucleon\n";
        dy << "# projectile   = " << beam_name << "\n";
        dy << "# current      = " << current_name << "\n";
        dy << "# units_sigma  = cm^2\n";
        dy << "# units_E      = GeV\n";
        dy << "# M_Z_GeV      = " << MZ << "\n";
        dy << "# dipole_model = " << model_name << "\n";
        dy << "# generated_by = sigma_nuN (PHASIS/dipole) commit "
           << PHASIS_GIT_HASH << "\n";
        dy << "# M_W_GeV      = " << MW << "\n";
        dy << "# xF3_from     = " << (useF3 ? pdf_set : std::string("nenhum")) << "\n";
        dy << "# Q2min_GeV2   = " << Q2min << "\n";
        dy << "# y_min_grade  = " << y_lo
           << "   (= Q2min/(s x_max) na energia mais alta da tabela)\n";
        dy << "# nE = " << NE << "   nY = " << nY << "\n";
        dy << "#\n";
        dy << "# CORTE CINEMATICO: dsigma/dy = 0 para y < Q2min/(s x_max),\n";
        dy << "# que cresce quando E cai. Os zeros nas energias baixas sao\n";
        dy << "# consequencia de Q^2 >= Q2min, nao lacunas da tabela.\n";
        dy << "#\n";
        dy << "# fracao de sigma recuperada por INT dy (dsigma/dy) sobre esta grade:\n";

        std::vector<double> Es(NE), ys(nY);
        for (int i = 0; i < NE; ++i) {
            const double logE = (NE == 1) ? logEmin
                : logEmin + i*(logEmax - logEmin)/(NE - 1);
            Es[i] = std::pow(10.0, logE);
        }
        for (int b = 0; b < nY; ++b) {
            ys[b] = y_lo*std::pow(1.0/y_lo, static_cast<double>(b)/(nY - 1));
        }

        std::vector<std::vector<double>> D(NE, std::vector<double>(nY, 0.0));

        // Paralelo SO sem LHAPDF. LHAPDF::GridPDF guarda cache mutavel
        // dentro do objeto, entao compartilhar um PDF entre threads e
        // corrida de dados -- e o sintoma seria ruido pequeno e
        // irreprodutivel na tabela, exatamente o tipo de coisa que
        // levaria dias para diagnosticar depois.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) if(!useF3)
#endif
        for (int i = 0; i < NE; ++i) {
            for (int b = 0; b < nY; ++b) {
                D[i][b] = dsigma_dy(Es[i], ys[b], current, table, nodesX,
                                    useF3, weakSF.get(), beam, Q2min)
                          *GeVminus2_to_cm2;
            }
        }

        for (int i = 0; i < NE; ++i) {
            // trapezio em y linear, que e o que o leitor do PHASIS usa
            double acc = 0.0;
            for (int b = 0; b + 1 < nY; ++b) {
                acc += 0.5*(D[i][b] + D[i][b+1])*(ys[b+1] - ys[b]);
            }
            int nq = 0;
            const double sig = sigmaNuN(Es[i], current, table, nodesQ, nodesX,
                                        useF3, weakSF.get(), beam, Q2min, &nq)
                               *GeVminus2_to_cm2;
            if (i % std::max(1, NE/12) == 0 || i == NE-1) {
                dy << "#   E = " << Es[i] << " GeV : " << (sig > 0.0 ? acc/sig : 0.0)
                   << "   (y_min cinematico = " << yMinKinematic(Es[i], Q2min) << ")\n";
            }
        }

        dy << "# E_GeV y dsigma_dy_cm2\n";
        for (int i = 0; i < NE; ++i) {
            for (int b = 0; b < nY; ++b) {
                dy << Es[i] << " " << ys[b] << " " << D[i][b] << "\n";
            }
        }
        std::cout << "Arquivo gerado: " << dy_file << "\n";
    }

    if (table.clampedQueries() > 0) {
        std::cerr << "AVISO: " << table.clampedQueries()
                  << " consultas fora da faixa da tabela foram truncadas.\n";
    }

    std::cout << "Arquivo gerado: " << output_file << "\n";
    return 0;
}
