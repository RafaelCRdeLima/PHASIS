#include "sigma_nuN_core.hpp"

#include <cmath>

namespace dipole {

double largeXFactor(double x)
{
    if (x >= 1.0) return 0.0;
    if (x <= 0.0) return 0.0;
    return std::pow(1.0 - x, 7.0);
}

// =====================================================
// F4 - numero de nos proporcional a faixa de integracao
//
// O codigo antigo usava NlogQ = 16 FIXO para uma faixa
// ln(1) ate ln(0.999*2*M_N*E), que cresce com a energia. Em 1e6 GeV
// isso da ~5 nos cobrindo o pico em Q^2 ~ M_W^2; em 1e14 GeV, ~2.
// O resultado parava de crescer e ficava ruidoso, e uma perturbacao
// de 1e-7 em E movia sigma em 2%.
//
// Fixando a DENSIDADE de nos por decada, o espacamento em ln Q^2 nao
// muda com a energia e o pico e sempre resolvido igual.
// =====================================================
int nodesForRange(double log_lo, double log_hi,
                         double per_decade, int n_min)
{
    const double decades = (log_hi - log_lo)/std::log(10.0);

    int n = static_cast<int>(std::ceil(per_decade*decades));
    if (n < n_min) n = n_min;
    if (n % 2 != 0) ++n;

    return n;
}

namespace {

// Massa do boson trocado. E a UNICA constante que distingue as duas
// correntes na formula mestra -- todo o resto da diferenca ja esta
// dentro de F_2, F_L e xF_3.
double bosonMass2(CurrentType current)
{
    return (current == CurrentType::NC) ? MZ*MZ : MW*MW;
}

// Nucleo comum. Recebe x, y e Q^2 JA CONSISTENTES entre si; quem chama
// decide qual dos dois deriva do outro.
double d2sigma_core(
    double Enu, double x, double y, double Q2,
    CurrentType current,
    const StructureTable& table,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam
)
{
    // A tabela guarda F_T e F_L crus; o (1-x)^7 (regra de contagem de
    // constituintes, Kutak-Kwiecinski) e aplicado aqui para manter a
    // tabela suave.
    const StructureTL raw = table.at(x, Q2);
    const double lx = largeXFactor(x);

    const double FT = lx*raw.FT;
    const double FL = lx*raw.FL;

    const double MV2 = bosonMass2(current);
    const double propagator = std::pow(MV2/(Q2 + MV2), 2);
    const double prefactor  = GF*GF*MN*Enu/pi * propagator;

    double xF3 = 0.0;
    if (useF3 && weakSF != nullptr) {
        xF3 = (current == CurrentType::NC)
            ? weakSF->xF3_NC_isoscalar(x, Q2, beam)
            : weakSF->xF3_CC_isoscalar(x, Q2, beam);
    }

    // O sinal de nu/nubar ja esta dentro de xF3_NC_isoscalar; para o CC
    // ele e aplicado aqui, como sempre foi. Manter as duas convencoes
    // seria pedir um erro de sinal, entao o CC continua com o fator
    // explicito e o NC devolve o valor ja assinado.
    const double f3_sign = (current == CurrentType::NC)
        ? 1.0
        : ((beam == weak::BeamType::Neutrino) ? 1.0 : -1.0);

    const double bracket =
        0.5*(1.0 + std::pow(1.0-y, 2))*(FT + FL)
        - 0.5*y*y*FL
        + f3_sign*y*(1.0 - 0.5*y)*xF3;

    return prefactor*bracket;
}

// y = 1 e o EXTREMO CINEMATICO (lepton de saida com energia nula), nao
// um ponto proibido: o integrando ali vale
//
//     prefator * [ F2/2 - F_L/2 + xF3/2 ]
//
// que e finito e da ordem do resto. Zera-lo -- como fazia o teste
// `y >= 1.0` -- punha um DEGRAU exatamente no extremo inferior da
// quadratura em x, porque o laco integra a partir de x = Q^2/s, que e
// precisamente y = 1. Simpson com o valor do no de borda errado erra em
// O(h), nao em O(h^4): o erro deixa de cair como N^-4 e passa a cair
// como N^-1.
//
// Foi o que se mediu. Triplicar a densidade de nos reduzia o ruido de
// sigma por 3.2, nao por 81; e aumentar a tabela de F de 121 para 241
// nos nao mudava NADA (as duas corridas sairam identicas), provando que
// a fonte estava na quadratura, nao na interpolacao. Depois do conserto
// o ruido caiu por um fator 670 e as duas densidades de nos concordam
// em 5e-7.
//
// Mesma classe dos dois bugs de borda do PHASIS: fronteira dura
// avaliada exatamente num extremo de integracao.
bool clampY(double& y)
{
    if (!(y > 0.0)) return false;
    if (y > 1.0) {
        // Fora da cinematica de verdade, ou so o arredondamento de
        // exp(log(Q2/s)) no no de borda?
        if (y > 1.0 + 1.0e-9) return false;
        y = 1.0;
    }
    return true;
}

} // namespace

// ---------------------------------------------------------------------

double d2sigma_dxdQ2_param(
    double Enu, double x, double Q2,
    CurrentType current,
    const StructureTable& table,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam
)
{
    if (x <= 0.0 || x >= 1.0) return 0.0;

    const double s = 2.0*MN*Enu;
    double y = Q2/(x*s);
    if (!clampY(y)) return 0.0;

    return d2sigma_core(Enu, x, y, Q2, current, table, useF3, weakSF, beam);
}

double d2sigma_dxdy_param(
    double Enu, double x, double y,
    CurrentType current,
    const StructureTable& table,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam
)
{
    if (x <= 0.0 || x >= 1.0) return 0.0;
    if (!(y > 0.0) || y > 1.0) return 0.0;

    const double s  = 2.0*MN*Enu;
    const double Q2 = x*y*s;
    if (!(Q2 > 0.0)) return 0.0;

    return d2sigma_core(Enu, x, y, Q2, current, table, useF3, weakSF, beam);
}

double d2sigma_dxdy_CC(
    double Enu, double x, double Q2,
    const StructureTable& table,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam
)
{
    return d2sigma_dxdQ2_param(Enu, x, Q2, CurrentType::CC,
                               table, useF3, weakSF, beam);
}

// =====================================================================

double sigmaNuN(
    double Enu,
    CurrentType current,
    const StructureTable& table,
    double nodesPerDecadeQ,
    double nodesPerDecadeX,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam,
    double Q2min,
    int* nodesQout
)
{
    const double s = 2.0*MN*Enu;
    const double Q2max = 0.999*s;

    if (Q2max <= Q2min) return 0.0;

    const double logQmin = std::log(Q2min);
    const double logQmax = std::log(Q2max);

    const int NlogQ = nodesForRange(logQmin, logQmax, nodesPerDecadeQ, 16);
    if (nodesQout) *nodesQout = NlogQ;

    const double xmax = 0.999999;

    // O numero de nos em x tem de ser FIXO ao longo do laco de Q^2.
    //
    // Se for recalculado por no (a faixa e [Q^2/s, xmax], que encolhe
    // conforme Q^2 cresce), ele salta de 2 em 2 ao longo do laco, e cada
    // salto e uma descontinuidade no integrando EXTERNO: o erro de
    // convergencia da integral interna muda de degrau. Deslocar os nos
    // externos por 1e-7 fazia alguns cruzarem esses degraus, e sigma
    // pulava 0.2-0.6% - o que sobrava do ruido diagnosticado na secao 0.
    //
    // Fixando pela faixa mais larga (a de Q^2 = Q2min), o integrando
    // externo volta a ser suave e as faixas estreitas ficam apenas
    // super-resolvidas, o que so custa tempo.
    const int Nlogx = nodesForRange(std::log(Q2min/s), std::log(xmax),
                                    nodesPerDecadeX, 16);

    auto integrandLogQ = [&](double logQ2)
    {
        const double Q2 = std::exp(logQ2);

        const double xmin = Q2/s;

        if (xmin >= xmax) return 0.0;

        const double logxmin = std::log(xmin);
        const double logxmax = std::log(xmax);

        auto integrandLogx = [&](double logx)
        {
            const double x = std::exp(logx);

            return d2sigma_dxdQ2_param(
                       Enu, x, Q2, current, table, useF3, weakSF, beam
                   )/s;
        };

        return Q2*simpson(integrandLogx, logxmin, logxmax, Nlogx);
    };

    return simpson(integrandLogQ, logQmin, logQmax, NlogQ);
}

double sigmaNuN_CC(
    double Enu,
    const StructureTable& table,
    double nodesPerDecadeQ,
    double nodesPerDecadeX,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam,
    double Q2min,
    int* nodesQout
)
{
    return sigmaNuN(Enu, CurrentType::CC, table, nodesPerDecadeQ,
                    nodesPerDecadeX, useF3, weakSF, beam, Q2min, nodesQout);
}

// =====================================================================

double yMinKinematic(double Enu, double Q2min, double xmax)
{
    const double s = 2.0*MN*Enu;
    return Q2min/(s*xmax);
}

double dsigma_dy(
    double Enu,
    double y,
    CurrentType current,
    const StructureTable& table,
    double nodesPerDecadeX,
    bool useF3,
    const weak::WeakStructureFunctions* weakSF,
    weak::BeamType beam,
    double Q2min
)
{
    if (!(y > 0.0) || y > 1.0) return 0.0;

    const double s = 2.0*MN*Enu;
    const double xmax = 0.999999;

    // Q^2 = x y s >= Q2min  =>  x >= Q2min/(y s)
    const double xmin = Q2min/(y*s);
    if (xmin >= xmax) return 0.0;

    // Mesma regra de densidade de nos por decada da rota (Q^2, x). Aqui
    // a faixa em x encolhe quando y cresce, e o numero de nos acompanha:
    // nao ha laco externo para tornar o salto de 2 em 2 visivel, entao a
    // razao para congelar N nao se aplica.
    const int Nlogx = nodesForRange(std::log(xmin), std::log(xmax),
                                    nodesPerDecadeX, 16);

    auto integrandLogx = [&](double logx)
    {
        const double x = std::exp(logx);
        // jacobiano de dx = x d(ln x)
        return x*d2sigma_dxdy_param(Enu, x, y, current, table,
                                    useF3, weakSF, beam);
    };

    return simpson(integrandLogx, std::log(xmin), std::log(xmax), Nlogx);
}

} // namespace dipole
