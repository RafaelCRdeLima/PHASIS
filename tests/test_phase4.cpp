// Fase 4, etapas 2 e 3: T20 (limite sem regeneracao), T21 (conservacao),
// T22 (solucao analitica com perda e ganho juntos).
// =====================================================================
// Fase 4 -- cascata NC. Etapas 1 a 3: identidade discreta, T20, T21, T22.
// (T23-T26 pendentes.)
// =====================================================================

#include "phasis/cascade.hpp"
#include "phasis/cross_section.hpp"
#include "phasis/metric.hpp"
#include "phasis/trace.hpp"
#include "phasis/units.hpp"
#include "phasis/density.hpp"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <vector>
using namespace phasis;

static int falhas=0, total=0;
static void chk(const char* n,double g,double r,double t){
    ++total; const double d=(std::fabs(r)>0?std::fabs(r):1.0);
    const double rel=std::fabs(g-r)/d; const bool ok=rel<=t&&std::isfinite(g);
    if(!ok)++falhas;
    std::printf("    %-44s %-18.11g ref %-18.11g rel %8.2e  %s\n",n,g,r,rel,ok?"OK":"<-- FALHA");
}

int main(){
    const Minkowski flat;
    const double R = 1.0e10, rho0 = 1.0e-3;
    const UniformBall laje(rho0, R);
    const double L = 2.0*R;                       // travessia diametral
    const double nN = units::nucleons_per_gram*rho0;

    std::printf("\n=====================================================\n");
    std::printf(" PHASIS -- Fase 4: cascata NC (etapas 2 e 3)\n");
    std::printf("=====================================================\n");

    // ---------------------------------------------------------------
    std::printf("\n  ETAPA 1  identidade discreta  soma_i G_ij + G_leak_j == 1\n\n");
    std::printf("  %-14s %-10s %-8s %-16s %s\n","kernel","bins/dec","N","pior residuo","leak(topo)");
    double pior_id = 0.0;

    
    for (int npd : {5, 10, 20, 40}) {
        const int N = npd*8;                      // 1e3 a 1e11 GeV
        EnergyGrid g(1.0e3, 1.0e11, N);

        {   ConstantY xs(1.0e-33, 2.4e-33, 1.0e6, 0.36, 0.36);
            CascadeKernel k(g, xs, 1.0e-12);
            std::printf("  %-14s %-10d %-8d %-16.3e %.6f\n","ConstantY",npd,N,
                        k.worst_residual(), k.G_leak(N-1));
            pior_id = std::max(pior_id, k.worst_residual());
        }
        {   PowerLawY xs(1.0e-33, 2.4e-33, 1.0e6, 0.36, 0.36, 1.0, 1.0e-4);
            CascadeKernel k(g, xs, 1.0e-12);
            std::printf("  %-14s %-10d %-8d %-16.3e %.6f\n","PowerLawY b=1",npd,N,
                        k.worst_residual(), k.G_leak(N-1));
            pior_id = std::max(pior_id, k.worst_residual());
        }
    }
    std::printf("\n  Z_discrete vs Z analitico (ConstantY, Z = 1/gamma)\n\n");
    std::printf("  %-10s %-10s %-16s %-16s %s\n","gamma","bins/dec","Z_disc","Z_analitico","desvio rel");
    for (double gam : {1.5, 2.0, 2.7}) {
        for (int npd : {10, 40, 160}) {
            EnergyGrid g(1.0e3, 1.0e11, npd*8);
            ConstantY xs(1.0e-33, 2.4e-33, 1.0e6, 0.36, 0.36);
            CascadeKernel k(g, xs, 1.0e-12);
            const double zd = k.Z_discrete(gam), za = ConstantY::Z_analytic(gam);
            std::printf("  %-10.1f %-10d %-16.10f %-16.10f %.3e\n",
                        gam,npd,zd,za,std::fabs(zd-za)/za);
        }
    }
    chk("identidade: pior residuo", pior_id, 0.0, 1.0e-12);

    std::printf("\n  T20  limite sem regeneracao: dsigma_NC/dy = 0\n");
    {
        EnergyGrid g(1.0e5, 1.0e11, 60);
        ZeroNC xs(2.4e-33, 1.0e6, 0.0);           // sigma_cc constante
        CascadeKernel k(g, xs, 1.0e-12);

        std::vector<double> phi0(g.n(), 0.0);
        // um bin so, para comparar com exp(-tau) puro
        const int jb = 40;
        phi0[jb] = 1.0;

        Ray ray; ray.E_inf_GeV = g.center(jb); ray.b_cm = 0.0;
        IntegratorOpts o; o.ode_rel_tol = 1.0e-13;
        const CascadeResult r = transport_cascade(phi0, k, ray, flat, laje, o);

        const double tau = nN*xs.sigma_cc(g.center(jb))*L;
        chk("phi/phi0 == exp(-tau)", r.phi[jb], std::exp(-tau), 1.0e-11);
        chk("nenhum bin abaixo recebeu fluxo", r.phi[jb-1], 0.0, 0.0);
        chk("vazamento nulo", r.leakage_total, 0.0, 0.0);
        chk("coluna", r.column_g_cm2, rho0*L, 1.0e-10);
    }

    // ---------------------------------------------------------------
    std::printf("\n  T21  conservacao de numero: sigma_CC = 0, so NC\n");
    {
        EnergyGrid g(1.0e5, 1.0e11, 90);
        ConstantY xs(1.0e-33, 0.0, 1.0e6, 0.0, 0.0);   // sigma_cc = 0
        CascadeKernel k(g, xs, 1.0e-12);

        // espectro inicial em lei de potencia
        std::vector<double> phi0(g.n());
        for (int i=0;i<g.n();++i) phi0[i]=g.powerlaw_bin_integral(i,2.0);
        const double N0 = std::accumulate(phi0.begin(), phi0.end(), 0.0);

        Ray ray; ray.E_inf_GeV=1.0e8; ray.b_cm=0.0;
        IntegratorOpts o; o.ode_rel_tol=1.0e-13;

        std::printf("    %-14s %-16s %-16s %s\n","fracao de L","soma phi","+ vazamento","desvio rel");
        double pior=0;
        for (double frac : {0.25, 0.5, 1.0}) {
            const UniformBall parcial(rho0, frac*R);
            const CascadeResult r = transport_cascade(phi0,k,ray,flat,parcial,o);
            const double Nf = std::accumulate(r.phi.begin(), r.phi.end(), 0.0);
            const double tot = Nf + r.leakage_total;
            const double dev = std::fabs(tot-N0)/N0;
            pior = std::max(pior, dev);
            std::printf("    %-14.2f %-16.10e %-16.10e %.3e\n", frac, Nf, tot, dev);
        }
        chk("conservacao em todos os pontos", pior, 0.0, 1.0e-12);
    }

    // ---------------------------------------------------------------
    std::printf("\n  T22  solucao analitica: perda E ganho juntos\n");
    {
        const double snc = 1.0e-33, scc = 5.0e-34;
        std::printf("    %-8s %-10s %-16s %-16s %-12s %s\n",
                    "gamma","bins/dec","phi(L)/phi(0)","exp(-n sig_ef L)","desvio","sig_ef/sig_ef,an");
        for (double gam : {1.5, 2.0, 2.7}) {
            for (int npd : {15, 60}) {
                EnergyGrid g(1.0e2, 1.0e14, npd*12);
                ConstantY xs(snc, scc, 1.0e6, 0.0, 0.0);   // sem dependencia em E
                CascadeKernel k(g, xs, 1.0e-12);

                std::vector<double> phi0(g.n());
                for(int i=0;i<g.n();++i) phi0[i]=g.powerlaw_bin_integral(i,gam);

                Ray ray; ray.E_inf_GeV=1.0e8; ray.b_cm=0.0;
                IntegratorOpts o; o.ode_rel_tol=1.0e-13;
                const CascadeResult r=transport_cascade(phi0,k,ray,flat,laje,o);

                // bin de referencia, longe das bordas
                const int ib = g.n()/2;
                const double razao = r.phi[ib]/phi0[ib];

                const double Zd  = k.Z_discrete(gam);
                const double sef = scc + snc*(1.0 - Zd);
                const double esp = std::exp(-nN*sef*L);

                const double sefa = scc + snc*(1.0 - ConstantY::Z_analytic(gam));

                std::printf("    %-8.1f %-10d %-16.10e %-16.10e %-12.3e %.8f\n",
                            gam,npd,razao,esp,std::fabs(razao-esp)/esp,sef/sefa);
                char nm[80]; std::snprintf(nm,sizeof(nm),"gamma=%.1f npd=%d",gam,npd);
                chk(nm, razao, esp, 1.0e-10);
            }
        }
    }

    std::printf("\n=====================================================\n");
    std::printf("  %d verificacoes, %d falhas\n", total, falhas);
    std::printf("  RESULTADO: %s\n", falhas==0?"OK":"FALHA");
    std::printf("=====================================================\n\n");
    return falhas==0?0:1;
}
