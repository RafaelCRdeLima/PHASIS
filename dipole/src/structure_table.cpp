#include "structure_table.hpp"

#include <chrono>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <stdexcept>

namespace dipole {

namespace {

// Sentinela para "F <= 0 neste no". ln(F) nunca chega perto disso.
constexpr double LOG_ZERO = -1.0e30;

double logOrSentinel(double v)
{
    return (v > 0.0 && std::isfinite(v)) ? std::log(v) : LOG_ZERO;
}

// Interpolacao bicubica (Catmull-Rom) em (ln x, ln Q^2) sobre ln F.
//
// Bilinear e C^0 mas nao C^1: a derivada salta em cada linha da grade.
// A integracao em sigma_nuN amostra a tabela com nos de Simpson que se
// deslocam quando a energia muda, e a cada deslocamento alguns nos
// cruzam essas quinas. Medido: residuo de 0.46% rms em sigma(E) so por
// causa disso.
//
// Catmull-Rom e C^1 e local (estencil 4x4), o que remove as quinas sem
// custo de montagem nem de memoria.
double catmullRom(double t, double p0, double p1, double p2, double p3)
{
    return 0.5*( 2.0*p1
               + (-p0 + p2)*t
               + (2.0*p0 - 5.0*p1 + 4.0*p2 - p3)*t*t
               + (-p0 + 3.0*p1 - 3.0*p2 + p3)*t*t*t );
}

} // namespace

StructureTable::StructureTable(
    DipoleModelId model,
    const std::vector<Channel>& channels,
    const GBWParameters& gbw,
    const IIMParameters& iim,
    const QuadratureGrid& quad,
    const TableSpec& spec,
    const QuarkMasses& masses,
    bool verbose
)
    : model_(model), channels_(channels), gbw_(gbw), iim_(iim),
      quad_(quad), spec_(spec), masses_(masses)
{
    if (spec_.Nx < 2 || spec_.NQ < 2) {
        throw std::runtime_error("TableSpec precisa de ao menos 2 pontos por eixo.");
    }
    if (!(spec_.xMin > 0.0) || !(spec_.xMax > spec_.xMin)) {
        throw std::runtime_error("Faixa de x invalida em TableSpec.");
    }
    if (!(spec_.Q2Min > 0.0) || !(spec_.Q2Max > spec_.Q2Min)) {
        throw std::runtime_error("Faixa de Q2 invalida em TableSpec.");
    }

    const auto t0 = std::chrono::steady_clock::now();

    lx_.resize(spec_.Nx);
    lq_.resize(spec_.NQ);

    for (int i = 0; i < spec_.Nx; ++i) {
        lx_[i] = std::log(spec_.xMin)
               + i*(std::log(spec_.xMax) - std::log(spec_.xMin))/(spec_.Nx - 1);
    }
    for (int j = 0; j < spec_.NQ; ++j) {
        lq_[j] = std::log(spec_.Q2Min)
               + j*(std::log(spec_.Q2Max) - std::log(spec_.Q2Min))/(spec_.NQ - 1);
    }

    FT_.assign(static_cast<std::size_t>(spec_.Nx)*spec_.NQ, LOG_ZERO);
    FL_.assign(static_cast<std::size_t>(spec_.Nx)*spec_.NQ, LOG_ZERO);

    if (verbose) {
        std::fprintf(stderr, "  tabelando F_T e F_L em %d x %d pontos...\n",
                     spec_.Nx, spec_.NQ);
    }

    int done = 0;
    int last_decile = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < spec_.Nx; ++i) {

        const double x = std::exp(lx_[i]);

        for (int j = 0; j < spec_.NQ; ++j) {

            const double Q2 = std::exp(lq_[j]);
            const StructureTL F = computeAt(x, Q2);

            const std::size_t k = static_cast<std::size_t>(i)*spec_.NQ + j;
            FT_[k] = logOrSentinel(F.FT);
            FL_[k] = logOrSentinel(F.FL);
        }

        if (verbose) {
#ifdef _OPENMP
#pragma omp atomic
#endif
            ++done;

            const int pct = static_cast<int>(100.0*done/spec_.Nx);
            if (pct/10 > last_decile) {
                last_decile = pct/10;
                std::fprintf(stderr, "    %3d%%\n", 10*last_decile);
            }
        }
    }


    const auto t1 = std::chrono::steady_clock::now();
    build_seconds_ =
        std::chrono::duration<double>(t1 - t0).count();
}

StructureTL StructureTable::computeAt(double x, double Q2) const
{
    StructureTL total;

    for (const Channel& ch : channels_) {

        Parameters wf = makeCCParameters(ch.first, ch.second, masses_);

        const StructureTL F = (model_ == DipoleModelId::GBW)
            ? structureGBW(x, Q2, wf, gbw_, quad_)
            : structureIIM(x, Q2, wf, iim_, quad_);

        // Convencao DIS: remove alphaEW global das wavefunctions.
        // Ver CAMPANHA_CORRECAO.md secao 4, Q1.
        total.FT += F.FT/wf.alphaEW;
        total.FL += F.FL/wf.alphaEW;
    }

    return total;
}

StructureTL StructureTable::exact(double x, double Q2) const
{
    return computeAt(x, Q2);
}

StructureTL StructureTable::at(double x, double Q2) const
{
    if (!(x > 0.0) || !(Q2 > 0.0)) return StructureTL{};

    double lx = std::log(x);
    double lq = std::log(Q2);

    const double lx0 = lx_.front(), lx1 = lx_.back();
    const double lq0 = lq_.front(), lq1 = lq_.back();

    if (lx < lx0 || lx > lx1 || lq < lq0 || lq > lq1) {
        ++clamped_;
        if (lx < lx0) lx = lx0;
        if (lx > lx1) lx = lx1;
        if (lq < lq0) lq = lq0;
        if (lq > lq1) lq = lq1;
    }

    const double hx = (lx1 - lx0)/(spec_.Nx - 1);
    const double hq = (lq1 - lq0)/(spec_.NQ - 1);

    int i = static_cast<int>((lx - lx0)/hx);
    int j = static_cast<int>((lq - lq0)/hq);

    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i > spec_.Nx - 2) i = spec_.Nx - 2;
    if (j > spec_.NQ - 2) j = spec_.NQ - 2;

    const double t = (lx - lx_[i])/hx;
    const double u = (lq - lq_[j])/hq;

    // estencil 4x4 em torno de (i,j), com indices grampeados na borda
    auto amostra = [&](const std::vector<double>& F, int di, int dj) {
        int a = i + di, b = j + dj;
        if (a < 0) a = 0;
        if (b < 0) b = 0;
        if (a > spec_.Nx - 1) a = spec_.Nx - 1;
        if (b > spec_.NQ - 1) b = spec_.NQ - 1;
        return F[static_cast<std::size_t>(a)*spec_.NQ + b];
    };

    auto interpola = [&](const std::vector<double>& F) {
        double col[4];
        for (int di = -1; di <= 2; ++di) {
            const double q0 = amostra(F, di, -1);
            const double q1 = amostra(F, di,  0);
            const double q2 = amostra(F, di,  1);
            const double q3 = amostra(F, di,  2);

            // sentinela em qualquer ponto do estencil -> nao ha valor
            if (q0 <= LOG_ZERO/2 || q1 <= LOG_ZERO/2
             || q2 <= LOG_ZERO/2 || q3 <= LOG_ZERO/2) return 0.0;

            col[di + 1] = catmullRom(u, q0, q1, q2, q3);
        }
        return std::exp(catmullRom(t, col[0], col[1], col[2], col[3]));
    };

    StructureTL out;
    out.FT = interpola(FT_);
    out.FL = interpola(FL_);
    return out;
}

std::string StructureTable::describe() const
{
    std::ostringstream os;
    os << (model_ == DipoleModelId::GBW ? "GBW" : "IIM")
       << "  Nx=" << spec_.Nx << " NQ=" << spec_.NQ
       << "  x in [" << spec_.xMin << ", " << spec_.xMax << "]"
       << "  Q2 in [" << spec_.Q2Min << ", " << spec_.Q2Max << "]"
       << "  quad Nr=" << quad_.Nr << " Nz=" << quad_.Nz
       << "  canais=" << channels_.size();
    return os.str();
}

} // namespace dipole
