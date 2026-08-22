// =====================================================================
// Fase 3 -- testes de aceitacao T12..T18.
// (T19, nao-regressao, sao as suites das Fases 1 e 2 rodadas sem mudanca.)
// =====================================================================

#include "phasis/cross_section.hpp"
#include "phasis/density.hpp"
#include "phasis/metric.hpp"
#include "phasis/sweep.hpp"
#include "phasis/trace.hpp"
#include "phasis/units.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace phasis;

namespace {

constexpr double kPi = 3.141592653589793;

int falhas = 0, total = 0;

void check(const char* nome, double got, double ref, double tol_rel)
{
    ++total;
    const double d = (std::fabs(ref) > 0.0) ? std::fabs(ref) : 1.0;
    const double rel = std::fabs(got - ref)/d;
    const bool ok = (rel <= tol_rel) && std::isfinite(got);
    if (!ok) ++falhas;
    std::printf("    %-46s %-18.11g ref %-18.11g rel %8.2e  %s\n",
                nome, got, ref, rel, ok ? "OK" : "<-- FALHA");
}

void check_bool(const char* nome, bool got, bool ref)
{
    ++total;
    const bool ok = (got == ref);
    if (!ok) ++falhas;
    std::printf("    %-46s %-18s ref %-18s                %s\n",
                nome, got ? "true" : "false", ref ? "true" : "false",
                ok ? "OK" : "<-- FALHA");
}

IntegratorOpts opt(bool ode = false, double rel = 1.0e-13)
{
    IntegratorOpts o;
    o.rel_tol     = rel;
    o.ode_rel_tol = rel;
    o.max_depth   = 45;
    o.force_ode   = ode;
    return o;
}

// Limite H -> infinito do FlaredThinDisk: rho = rho0 (r sin(theta)/R0)^(-p).
// NAO e PowerLawHalo, a nao ser no plano do disco -- ver T16.
class FlatDiskLimit final : public DensityProfile {
public:
    FlatDiskLimit(double rho0, double R0, double p, double r_in, double r_out)
        : rho0_(rho0), R0_(R0), p_(p), r_in_(r_in), r_out_(r_out) {}
    double rho(double r, double th) const override {
        if (r < r_in_ || r > r_out_) return 0.0;
        const double R = r*std::sin(th);
        if (!(R > 0.0)) return 0.0;
        return rho0_*std::pow(R/R0_, -p_);
    }
    double r_support_min() const override { return r_in_; }
    double r_support_max() const override { return r_out_; }
    bool   is_spherical()  const override { return false; }
    std::string name() const override { return "FlatDiskLimit"; }
private:
    double rho0_, R0_, p_, r_in_, r_out_;
};

} // namespace

int main()
{
    std::printf("\n=====================================================\n");
    std::printf(" PHASIS -- Fase 3: dependencia em theta e paralelismo\n");
    std::printf("=====================================================\n");

    const double r_s = 1.0e6;
    const Schwarzschild schw(r_s);
    const Minkowski flat;

    // =================================================================
    std::printf("\n  T12  equivalencia: rota de EDO vs quadratura das Fases 1-2\n");
    {
        const PowerLawHalo halo(1.0e3, 1.0e8, 2.0, 3.0*r_s, 1.0e4*r_s);
        const PowerLawCrossSection xs(1.0e-33, 1.0e6, 0.4);

        for (double frac : {3.0, 5.0, 20.0, 1.0e3}) {
            Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = frac*r_s;
            const double a = trace_ray(ray, schw, halo, xs, opt(false)).tau;
            const double b = trace_ray(ray, schw, halo, xs, opt(true)).tau;
            char nome[80];
            std::snprintf(nome, sizeof(nome), "tau  (b/r_s = %g)", frac);
            check(nome, b, a, 1.0e-12);
        }

        // caso quase-critico: e onde as duas rotas mais podem divergir
        const PowerLawHalo perto(1.0e3, 1.0e8, 2.0, 1.4*r_s, 1.0e4*r_s);
        Ray rc; rc.E_inf_GeV = 1.0e9;
        rc.b_cm = schw.b_crit()*(1.0 + 1.0e-8);
        const double a = trace_ray(rc, schw, perto, xs, opt(false, 1.0e-11)).tau;
        const double b = trace_ray(rc, schw, perto, xs, opt(true,  1.0e-11)).tau;
        std::printf("    near_critical: quadratura %.12g  EDO %.12g\n", a, b);
        check("tau  (near_critical)", b, a, 1.0e-9);
    }

    // =================================================================
    std::printf("\n  T13  identidade  dpsi/dl = b sqrt(f) / r^2\n");
    {
        double pior = 0.0;
        int n = 0;
        for (double frac : {3.5, 40.0, 5.0e3}) {
            const double b = frac*r_s;
            const double r_t = schw.r_turning(b, b);
            for (int k = 1; k <= 50; ++k) {
                const double r = r_t*(1.0 + 0.05*k);
                const RayDerivatives d = ray_derivatives(schw, r, r_t, b);
                const double lhs = d.dpsi_ds/d.dl_ds;
                const double rhs = b*std::sqrt(schw.f(r))/(r*r);
                pior = std::max(pior, std::fabs(lhs - rhs)/rhs);
                ++n;
            }
        }
        std::printf("    %d pontos em 3 raios\n", n);
        check("pior desvio da identidade", pior, 0.0, 1.0e-13);
    }

    // =================================================================
    std::printf("\n  T14  simetria z -> -z:  tau(i, psi_t) = tau(i, psi_t + pi)\n");
    {
        const FlaredThinDisk disco(1.0e3, 1.0e8, 1.0e6, 3.0*r_s, 1.0e4*r_s);
        const PowerLawCrossSection xs(1.0e-33, 1.0e6, 0.4);

        for (double i : {0.1, 0.7, 1.4}) {
            Ray a; a.E_inf_GeV = 1.0e9; a.b_cm = 30.0*r_s;
            a.inclination_rad = i; a.psi_turn_rad = 0.37;
            Ray b = a; b.psi_turn_rad = 0.37 + kPi;

            const Result ra = trace_ray(a, schw, disco, xs, opt());
            const Result rb = trace_ray(b, schw, disco, xs, opt());

            char nome[80];
            std::snprintf(nome, sizeof(nome), "tau  (i = %.1f rad)", i);
            check(nome, rb.tau, ra.tau, 1.0e-13);
        }
    }

    // =================================================================
    std::printf("\n  T15  assimetria dos ramos (pega o bug de dobrar um ramo)\n");
    {
        const FlaredThinDisk disco(1.0e3, 1.0e8, 3.0e6, 3.0*r_s, 3.0e3*r_s);
        const PowerLawCrossSection xs(1.0e-33, 1.0e6, 0.4);

        std::printf("    %-10s %-10s %-14s %-14s %s\n",
                    "i [rad]", "psi_t", "X_entrada", "X_saida", "|dif| relativa");

        // Dois regimes de proposito:
        //  * i pequeno: o raio fica perto do plano e AMBOS os ramos tem
        //    coluna apreciavel -- e o caso que melhor pega um bug de
        //    dobrar ramo, porque a diferenca e moderada em vez de
        //    "um ramo vale zero".
        //  * i maior: um ramo cruza o plano e o outro se afasta; a
        //    assimetria vira ordens de grandeza.
        int quantos = 0;
        for (double i : {0.05, 0.08, 0.6, 1.0}) {
            for (double pt : {0.25, 0.55}) {
                Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = 20.0*r_s;
                ray.inclination_rad = i; ray.psi_turn_rad = pt;
                const Result r = trace_ray(ray, schw, disco, xs, opt());

                const double m = 0.5*(r.column_inbound + r.column_outbound);
                const double d = std::fabs(r.column_inbound - r.column_outbound)/m;
                std::printf("    %-10.2f %-10.2f %-14.6e %-14.6e %.4f\n",
                            i, pt, r.column_inbound, r.column_outbound, d);
                ++total;
                if (!(d > 0.01)) { ++falhas; std::printf("      <-- FALHA: ramos iguais\n"); }
                else ++quantos;
            }
        }
        std::printf("    %d de 8 configuracoes com assimetria > 1%%\n", quantos);
    }

    // =================================================================
    std::printf("\n  T16  limite de espessura infinita\n");
    {
        const double rho0 = 1.0e3, R0 = 1.0e8, p = 15.0/8.0;
        const double rin = 3.0*r_s, rout = 1.0e4*r_s;
        const PowerLawCrossSection xs(1.0e-33, 1.0e6, 0.4);

        // (a) No plano do disco (i = 0) vale z = 0 e R = r, entao
        //     FlaredThinDisk e IDENTICO a PowerLawHalo -- para QUALQUER
        //     H, nao so no limite. Checagem forte do mapeamento (R,z).
        {
            const FlaredThinDisk d1(rho0, R0, 0.01*R0, rin, rout, p);
            const PowerLawHalo   d2(rho0, R0, p, rin, rout);
            Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = 30.0*r_s;
            ray.inclination_rad = 0.0; ray.psi_turn_rad = 0.4;
            check("i = 0: disco fino == PowerLawHalo (qualquer H)",
                  trace_ray(ray, schw, d1, xs, opt()).tau,
                  trace_ray(ray, schw, d2, xs, opt()).tau, 1.0e-12);
        }

        // (b) Fora do plano o limite H -> inf NAO e PowerLawHalo: como
        //     rho depende de R = r sin(theta), o limite e
        //     rho0 (r sin(theta)/R0)^(-p), que so coincide com
        //     PowerLawHalo em theta = pi/2. Testamos contra o limite
        //     correto.
        {
            const FlaredThinDisk d1(rho0, R0, 1.0e6*R0, rin, rout, p);
            const FlatDiskLimit  d2(rho0, R0, p, rin, rout);
            Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = 30.0*r_s;
            ray.inclination_rad = 0.9; ray.psi_turn_rad = 0.4;
            check("i = 0.9: H0/R0 = 1e6 -> limite (r sin th)^-p",
                  trace_ray(ray, schw, d1, xs, opt()).tau,
                  trace_ray(ray, schw, d2, xs, opt()).tau, 1.0e-6);
        }
    }

    // =================================================================
    std::printf("\n  T17  cruzamento em campo plano: forma fechada\n");
    {
        const double rho0 = 2.0, R0 = 1.0e8;
        const double rin = 1.0e7, rout = 1.0e9;
        const double sig = 1.0e-33;
        const PowerLawCrossSection xs(sig, 1.0, 0.0);
        const double NA = units::N_A;

        // i = 0: theta = pi/2, z = 0, R = r  =>  rho = rho0 (r/R0)^-p.
        // b = 0, dl = dr, integral direta:
        //   X = 2 rho0 R0^p [ r^(1-p)/(1-p) ]_{rin}^{rout}
        {
            const double p = 15.0/8.0;
            const FlaredThinDisk disco(rho0, R0, 0.01*R0, rin, rout, p);
            Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = 0.0;
            ray.inclination_rad = 0.0;
            const Result r = trace_ray(ray, flat, disco, xs, opt());

            const double X = 2.0*rho0*std::pow(R0, p)
                           *(std::pow(rout, 1.0-p) - std::pow(rin, 1.0-p))/(1.0 - p);
            check("p = 15/8, b = 0: coluna", r.column_density, X, 1.0e-10);
            check("p = 15/8, b = 0: tau",    r.tau, NA*sig*X,   1.0e-10);
        }

        // b > 0 com p = 2:  X = 2 rho0 R0^2 (1/b)[acos(b/rout) - acos(b/r_lo)]
        {
            const double p = 2.0;
            const FlaredThinDisk disco(rho0, R0, 0.01*R0, rin, rout, p);
            for (double b : {2.0e7, 2.0e8}) {
                Ray ray; ray.E_inf_GeV = 1.0e9; ray.b_cm = b;
                ray.inclination_rad = 0.0;
                const Result r = trace_ray(ray, flat, disco, xs, opt());

                const double r_lo = std::max(b, rin);
                const double X = 2.0*rho0*R0*R0*(1.0/b)
                               *(std::acos(b/rout) - std::acos(b/r_lo));
                char nome[80];
                std::snprintf(nome, sizeof(nome), "p = 2, b = %.0e: coluna", b);
                check(nome, r.column_density, X, 1.0e-9);
            }
        }
    }

    // =================================================================
    std::printf("\n  T18  determinismo do paralelismo\n");
    {
        const FlaredThinDisk disco(1.0e2, 1.0e8, 3.0e6, 3.0*r_s, 3.0e3*r_s);
        const PowerLawCrossSection xs(1.0e-33, 1.0e6, 0.4);

        std::vector<Ray> raios;
        raios.reserve(10000);
        for (int k = 0; k < 10000; ++k) {
            const double t = static_cast<double>(k)/9999.0;
            Ray r;
            r.E_inf_GeV       = 1.0e6*std::pow(1.0e6, t);
            r.b_cm            = r_s*(2.0 + 60.0*t);      // cruza b_crit
            r.inclination_rad = 0.2 + 1.2*t;
            r.psi_turn_rad    = 0.1 + 3.0*t;
            raios.push_back(r);
        }

        IntegratorOpts o = opt(false, 1.0e-9);
        std::vector<Result> ref;

#ifdef _OPENMP
        omp_set_num_threads(1);
#endif
        const SweepSummary sum1 = sweep(raios, schw, disco, xs, o, ref);
        std::printf("%s", ("    " + sum1.to_string()).c_str());

        bool identico = true;
        for (int nt : {2, 4, 8}) {
#ifdef _OPENMP
            omp_set_num_threads(nt);
#endif
            std::vector<Result> rr;
            sweep(raios, schw, disco, xs, o, rr);
            for (std::size_t k = 0; k < rr.size(); ++k) {
                const bool mesmo =
                    (rr[k].tau == ref[k].tau || (std::isinf(rr[k].tau) && std::isinf(ref[k].tau)))
                    && rr[k].tau_inbound  == ref[k].tau_inbound
                    && rr[k].tau_outbound == ref[k].tau_outbound
                    && rr[k].captured     == ref[k].captured;
                if (!mesmo) { identico = false; break; }
            }
            char nome[64];
            std::snprintf(nome, sizeof(nome), "identico bit a bit com %d threads", nt);
            check_bool(nome, identico, true);
        }
#ifdef _OPENMP
        omp_set_num_threads(omp_get_max_threads());
#endif
    }

    // =================================================================
    std::printf("\n  T18b  excecao nao escapa da regiao paralela\n");
    {
        // TableCrossSection lanca fora da faixa. Forcamos energias fora
        // dela: sem o try/catch dentro do laco isto seria comportamento
        // indefinido.
        const PowerLawHalo halo(1.0, 1.0e8, 2.0, 3.0*r_s, 1.0e4*r_s);

        struct XsecQueLanca final : CrossSection {
            double sigma_tot(double E) const override {
                if (E > 1.0e10) throw std::out_of_range("energia fora da faixa (teste)");
                return 1.0e-33;
            }
            std::string name() const override { return "XsecQueLanca"; }
        } xs;

        std::vector<Ray> raios;
        for (int k = 0; k < 2000; ++k) {
            Ray r; r.b_cm = 10.0*r_s;
            r.E_inf_GeV = (k % 2 == 0) ? 1.0e9 : 1.0e11;   // metade estoura
            raios.push_back(r);
        }

        std::vector<Result> res;
        const SweepSummary s = sweep(raios, schw, halo, xs, opt(false, 1.0e-8), res);
        std::printf("%s", ("    " + s.to_string()).c_str());
        check("raios com erro", static_cast<double>(s.n_errors), 1000.0, 0.0);
        check_bool("mensagem preservada",
                   !s.error_sample.empty()
                   && s.error_sample[0].find("fora da faixa") != std::string::npos, true);
    }

    std::printf("\n=====================================================\n");
    std::printf("  %d verificacoes, %d falhas\n", total, falhas);
    std::printf("  RESULTADO: %s\n", falhas == 0 ? "OK" : "FALHA");
    std::printf("=====================================================\n\n");
    return falhas == 0 ? 0 : 1;
}
