#ifndef DIPOLE_STRUCTURE_TABLE_HPP
#define DIPOLE_STRUCTURE_TABLE_HPP

#include <string>
#include <utility>
#include <vector>

#include "parameters.hpp"
#include "dipole_models.hpp"
#include "integrals.hpp"

namespace dipole {

enum class DipoleModelId { GBW, IIM };

// =====================================================
// Malha da tabela de F_T(x,Q^2) e F_L(x,Q^2)
//
// O ponto da tabela: F_{T,L} nao dependem da energia do neutrino,
// so de (x, Q^2). No codigo original a integral 2D em (r,z) rodava
// DENTRO dos lacos de x e Q^2, que por sua vez rodavam dentro do
// laco de energia -- refazendo o mesmo calculo 300 vezes.
//
// Tabelando uma vez e interpolando, o custo deixa de crescer com o
// numero de energias, e as grades de quadratura podem ser finas o
// suficiente para convergir (o que a F4 precisa).
// =====================================================
struct TableSpec {
    double xMin  = 1.0e-15;
    double xMax  = 1.0;
    double Q2Min = 1.0;
    double Q2Max = 1.0e15;

    int Nx = 121;   // pontos, nao intervalos
    int NQ = 121;
};

class StructureTable {
public:
    using Channel = std::pair<QuarkFlavor, QuarkFlavor>;

    StructureTable(
        DipoleModelId model,
        const std::vector<Channel>& channels,
        const GBWParameters& gbw,
        const IIMParameters& iim,
        const QuadratureGrid& quad,
        const TableSpec& spec,
        const QuarkMasses& masses,
        bool verbose,
        // CC: o dipolo junta DOIS sabores (ud~, cs~), e Channel e o par.
        // NC: o dipolo junta o mesmo sabor consigo (uu~, dd~, ...), e so
        // Channel::first e usado. Nao ha caminho separado para NC: a
        // mesma tabela, a mesma quadratura, so os acoplamentos e as
        // massas do par mudam.
        //
        // SEM VALOR PADRAO, de proposito. Ele tinha um, e foi assim que
        // a primeira tabela NC saiu errada: os canais eram uu, dd, ss,
        // cc mas os ACOPLAMENTOS continuaram CC (g_V=-1, g_A=1), porque
        // o argumento simplesmente nao foi passado. O resultado era
        // plausivel -- uma tabela suave, monotonica, sem excecao -- e
        // dava sigma_NC/sigma_CC = 2.49 em vez de 0.42. Um default
        // silencioso num parametro que muda a FISICA e um convite.
        CurrentType current
    );

    // F_T e F_L somados sobre os canais, ja divididos por alphaEW.
    // NAO aplica (1-x)^7: isso e responsabilidade de quem chama,
    // para manter a tabela suave.
    StructureTL at(double x, double Q2) const;

    // Mesma coisa, mas calculada direto, sem tabela. Para medir o
    // erro de interpolacao.
    StructureTL exact(double x, double Q2) const;

    CurrentType current() const { return current_; }
    const TableSpec& spec() const { return spec_; }
    double buildSeconds() const { return build_seconds_; }
    long clampedQueries() const { return clamped_; }
    std::string describe() const;

private:
    DipoleModelId model_;
    std::vector<Channel> channels_;
    CurrentType current_ = CurrentType::CC;
    GBWParameters gbw_;
    IIMParameters iim_;
    QuadratureGrid quad_;
    TableSpec spec_;
    QuarkMasses masses_;

    std::vector<double> lx_;   // ln x
    std::vector<double> lq_;   // ln Q^2
    std::vector<double> FT_;   // Nx * NQ
    std::vector<double> FL_;

    double build_seconds_ = 0.0;
    mutable long clamped_ = 0;

    StructureTL computeAt(double x, double Q2) const;
};

} // namespace dipole

#endif
