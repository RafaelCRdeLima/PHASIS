// Fase 4, etapas 2 e 3: T20 (limite sem regeneracao), T21 (conservacao),
// T22 (solucao analitica com perda e ganho juntos).
// =====================================================================
// Fase 4 -- cascata NC. Etapas 1 a 3: identidade discreta, T20, T21, T22.
// (T23-T26 pendentes.)
// =====================================================================

#include "phasis/cascade.hpp"
#include "phasis/cross_section.hpp"
#include "phasis/cross_section.hpp"
#include "phasis/metric.hpp"
#include <fstream>
#include <memory>
#include "phasis/trace.hpp"
#include "phasis/units.hpp"
#include "phasis/density.hpp"
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <vector>
using namespace phasis;

static constexpr double kPi = 3.141592653589793;
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


    // =================================================================
    std::printf("\n  BORDA  invariante generico, varrendo o registro de perfis\n\n");
    {
        // Duas ocorrencias do mesmo bug em fases diferentes: documentar
        // nao basta. O teste varre o REGISTRO, entao perfil novo entra
        // sozinho.
        //
        // O invariante que importa ao integrador nao e "rho > 0 na
        // borda" -- e que rho na borda seja o LIMITE PELO INTERIOR. Se
        // rho for nao-nula um fio para dentro e zero na borda, o estagio
        // de Runge-Kutta que cai ali ve um degrau que encolher o passo
        // nunca remove.
        const auto& reg = ProfileRegistry::instance().all();
        std::printf("    %d perfis registrados\n", static_cast<int>(reg.size()));

        const double thetas[] = { kPi/2, kPi/2 - 0.02, kPi/2 + 0.02, 0.6, 2.5 };

        for (const auto& e : reg) {
            const auto prof = e.second();
            const double a = prof->r_support_min();
            const double b = prof->r_support_max();

            bool ok = true;
            int testados = 0;
            for (double th : thetas) {
                for (int lado = 0; lado < 2; ++lado) {
                    const double borda = lado ? b : a;
                    const double eps   = 1.0e-9*std::max(borda, 1.0);
                    const double dentro = lado ? (borda - eps) : (borda + eps);

                    const double r_in  = prof->rho(dentro, th);
                    const double r_bd  = prof->rho(borda,  th);
                    if (r_in > 0.0) {
                        ++testados;
                        if (!(r_bd > 0.0)) ok = false;
                    }
                }
            }
            ++total; if (!ok) ++falhas;
            std::printf("    %-22s r em [%.3e, %.3e]  %d pontos  %s\n",
                        e.first.c_str(), a, b, testados, ok ? "OK" : "<-- FALHA");
        }
    }

    // =================================================================
    std::printf("\n  T23  exponencial de matriz, condicao inicial de BIN UNICO\n");
    {
        // Bin unico, nao lei de potencia: a lei de potencia e autovetor
        // EXATO do operador discreto, entao T22 testa uma direcao num
        // espaco de N dimensoes. Um erro que preserve sigma_eff para
        // leis de potencia mas erre a distribuicao entre bins passa
        // limpo por T22. Esta e a direcao que ele nao alcanca.
        EnergyGrid g(1.0e5, 1.0e11, 72);
        ConstantY xs(1.0e-33, 5.0e-34, 1.0e6, 0.0, 0.0);
        CascadeKernel k(g, xs, 1.0e-12);
        const int N = g.n();

        // M_ij = K_ij (j>i);  M_ii = K_ii - sigma_tot(E_i)
        std::vector<double> M(static_cast<std::size_t>(N)*N, 0.0);
        for (int j=0;j<N;++j){
            const double snc = xs.sigma_nc(g.center(j));
            for (int i=0;i<=j;++i) M[i*static_cast<std::size_t>(N)+j] = snc*k.G(i,j);
            M[j*static_cast<std::size_t>(N)+j] -= xs.sigma_cc(g.center(j)) + snc;
        }

        const double A = nN*L;
        for (double& v : M) v *= A;                // exponenciamos n L M

        // scaling-and-squaring com serie de Taylor
        int sq = 0; double nrm = 0.0;
        for (int i=0;i<N;++i){ double s=0; for(int j=0;j<N;++j) s+=std::fabs(M[i*static_cast<std::size_t>(N)+j]); nrm=std::max(nrm,s); }
        while (nrm > 0.5) { nrm *= 0.5; ++sq; }
        const double esc = std::pow(0.5, sq);
        for (double& v : M) v *= esc;

        auto mul=[&](const std::vector<double>& X,const std::vector<double>& Y){
            std::vector<double> Z(static_cast<std::size_t>(N)*N,0.0);
            for(int i=0;i<N;++i) for(int p=0;p<N;++p){
                const double x=X[i*static_cast<std::size_t>(N)+p]; if(x==0.0) continue;
                for(int j=0;j<N;++j) Z[i*static_cast<std::size_t>(N)+j]+=x*Y[p*static_cast<std::size_t>(N)+j];
            } return Z; };

        std::vector<double> E(static_cast<std::size_t>(N)*N,0.0), T(static_cast<std::size_t>(N)*N,0.0);
        for(int i=0;i<N;++i){ E[i*static_cast<std::size_t>(N)+i]=1.0; T[i*static_cast<std::size_t>(N)+i]=1.0; }
        for(int m=1;m<=24;++m){ T=mul(T,M); for(auto& v:T) v/=m; for(std::size_t q=0;q<E.size();++q) E[q]+=T[q]; }
        for(int s2=0;s2<sq;++s2) E=mul(E,E);

        IntegratorOpts o; o.ode_rel_tol=1.0e-13;
        Ray ray; ray.E_inf_GeV=1.0e8; ray.b_cm=0.0;

        std::printf("    %-8s %-14s %-14s %s\n","bin j","max |dif| rel","soma phi (EDO)","soma (matriz)");
        for (int jb : {N-2, N/2, 3}) {
            std::vector<double> phi0(N,0.0); phi0[jb]=1.0;
            const CascadeResult r=transport_cascade(phi0,k,ray,flat,laje,o);

            double pior=0.0, sE=0.0, sM=0.0;
            for(int i=0;i<N;++i){
                const double m=E[i*static_cast<std::size_t>(N)+jb];
                sE+=r.phi[i]; sM+=m;
                const double den=std::max(std::fabs(m),1.0e-30);
                if (std::fabs(m) > 1.0e-14) pior=std::max(pior,std::fabs(r.phi[i]-m)/den);
            }
            std::printf("    %-8d %-14.3e %-14.10f %.10f\n",jb,pior,sE,sM);
            char nm[64]; std::snprintf(nm,sizeof(nm),"bin unico j=%d vs exp(nLM)",jb);
            chk(nm,pior,0.0,1.0e-11);
        }
    }

    // =================================================================
    std::printf("\n  T23b  expansao de L curto: entradas INDIVIDUAIS de K_ij\n");
    {
        // A identidade de consistencia testa soma_i K_ij. Nada testava
        // cada K_ij. Com phi(0) = e_j e L -> 0,
        //     phi_i(L) = delta_ij + n L K_ij + O(L^2)
        // le cada entrada diretamente.
        EnergyGrid g(1.0e5, 1.0e11, 48);
        ConstantY xs(1.0e-33, 5.0e-34, 1.0e6, 0.0, 0.0);
        CascadeKernel k(g, xs, 1.0e-12);
        const int N=g.n();

        const double sig_max = xs.sigma_cc(g.center(N-1))+xs.sigma_nc(g.center(N-1));
        // n L sigma_max = 1e-10, nao 1e-6.
        //
        // A expansao e phi_i(L)/(nL) = K_ij + O(L), com erro RELATIVO da
        // ordem de n L sigma. Pedir 1e-8 de precisao com n L sigma = 1e-6
        // e incompativel: a truncagem sozinha ja vale 1e-6 (medido:
        // 3.4e-6). Com 1e-10 a truncagem fica em 1e-10 e a tolerancia de
        // 1e-8 passa a medir a implementacao, nao a expansao.
        const double Rfino = 1.0e-10/(2.0*nN*sig_max);
        const UniformBall fina(rho0, Rfino);
        const double Lf = 2.0*Rfino;

        IntegratorOpts o; o.ode_rel_tol=1.0e-14;
        Ray ray; ray.E_inf_GeV=1.0e8; ray.b_cm=0.0;

        std::printf("    n L sigma_max = %.2e  (truncagem O(L) da expansao)\n",
                    nN*Lf*sig_max);
        double pior=0.0; int checadas=0;
        for (int jb : {N-2, N/2, 5}) {
            std::vector<double> phi0(N,0.0); phi0[jb]=1.0;
            const CascadeResult r=transport_cascade(phi0,k,ray,flat,fina,o);
            const double snc=xs.sigma_nc(g.center(jb));
            for (int i=0;i<jb;++i){
                const double Kij = snc*k.G(i,jb);
                if (!(Kij > 1.0e-8*snc)) continue;         // ignora entradas nulas
                const double lido = r.phi[i]/(nN*Lf);
                pior=std::max(pior,std::fabs(lido-Kij)/Kij);
                ++checadas;
            }
        }
        std::printf("    %d entradas individuais lidas, pior desvio %.3e\n",checadas,pior);
        chk("cada K_ij lida pela expansao de L curto",pior,0.0,1.0e-8);
    }

    // =================================================================
    std::printf("\n  T24  convergencia em grade: ordem observada por JANELA\n");
    {
        // Previsao a checar: a grade em y induzida e y_k = 1 - rho^-k,
        // com espacamento ~ln(rho) perto de y = 0. Enquanto
        // ln(rho) > y_min o pico esta sub-resolvido e a ordem deve
        // degradar; quando ln(rho) < y_min o corte e resolvido e a ordem
        // deve VOLTAR para ~2. Espera-se uma TRANSICAO, nao um expoente
        // unico -- ajustar uma reta global daria um numero intermediario
        // sem significado.
        const double gam = 2.0;
        const int NP = 8;
        const int npds[NP] = {5,10,20,40,80,160,320,640};
        const double y_min = 0.05;

        for (int caso = 0; caso < 2; ++caso) {
            const bool picado = (caso == 1);
            std::printf("\n    %s\n", picado ? "PowerLawY  beta = 1,  y_min = 0.05"
                                             : "ConstantY  (suave)");
            std::printf("    %-8s %-11s %-15s %-12s %s\n",
                        "n/dec","ln(rho)","erro rel","ordem local","regime");

            double err[NP], lnr[NP];
            for (int q = 0; q < NP; ++q) {
                EnergyGrid g(1.0e6, 1.0e10, npds[q]*4);
                lnr[q] = std::log(g.ratio());
                double zd, za;
                if (picado) {
                    PowerLawY xs(1.0e-33,5.0e-34,1.0e6,0.0,0.0,1.0,y_min);
                    zd = z_discrete_row(g,xs,gam); za = xs.Z_analytic(gam);
                } else {
                    ConstantY xs(1.0e-33,5.0e-34,1.0e6,0.0,0.0);
                    zd = z_discrete_row(g,xs,gam); za = ConstantY::Z_analytic(gam);
                }
                err[q] = std::fabs(zd - za)/za;

                char ord[24] = "   -";
                if (q > 0) std::snprintf(ord,sizeof(ord),"%6.3f",
                                         std::log(err[q-1]/err[q])/std::log(2.0));
                const char* reg = !picado ? ""
                                : (lnr[q] < y_min ? "corte resolvido" : "sub-resolvido");
                std::printf("    %-8d %-11.4f %-15.4e %-12s %s\n",
                            npds[q], lnr[q], err[q], ord, reg);
            }

            if (!picado) {
                double m = 0.0;
                for (int q = 1; q < NP; ++q) m += std::log(err[q-1]/err[q])/std::log(2.0);
                chk("ConstantY: ordem media", m/(NP-1), 2.0, 0.05);
            } else {
                const double o_lo = std::log(err[0]/err[2])/std::log(4.0);
                const double o_hi = std::log(err[NP-3]/err[NP-1])/std::log(4.0);
                std::printf("\n      ordem sub-resolvida (n/dec  5 ->  20) : %.3f\n", o_lo);
                std::printf("      ordem resolvida     (n/dec 160-> 640) : %.3f\n", o_hi);
                std::printf("      transicao prevista em ln(rho) = y_min = %.3f  =>  n/dec = %.0f\n",
                            y_min, std::log(10.0)/y_min);
                ++total;
                if (!(o_hi > o_lo + 0.3)) {
                    ++falhas;
                    std::printf("      <-- FALHA: a ordem nao melhora ao resolver o corte\n");
                } else {
                    std::printf("      OK: ordem sobe de %.2f para %.2f ao cruzar y_min\n", o_lo, o_hi);
                }
            }
        }
    }

    // =================================================================
    std::printf("\n  T25  magnitude do redshift: M(E_loc) vs M congelada em E_inf\n");
    {
        // E_loc/E_inf = 1/sqrt(f(r_t)) <= sqrt(3) em Schwarzschild, porque
        // todo raio nao capturado tem r_t >= 1.5 r_s. Com sigma ~ E^0.36
        // a diferenca em P_surv tem de ficar em poucos por cento; uma
        // ordem de grandeza seria bug. O limite fisico serve de sanidade
        // numerica.
        //
        // A comparacao e com a MESMA geometria: freeze_redshift so troca
        // E_loc por E_inf no argumento de sigma. Comparar contra
        // Minkowski mudaria o caminho tambem.
        const double r_s = 1.0e6;
        const Schwarzschild schw(r_s);
        // Densidade escolhida para dar tau ~ 1. Com o halo fino de
        // antes (rho0 = 1e-3) a atenuacao total era 0.3%, e uma mudanca
        // de 11% em sigma aparecia como 0.03% no fluxo: o teste passava
        // sem medir nada. Aqui a diferenca fica na escala que o limite
        // sqrt(3) preve.
        const PowerLawHalo halo(3.0, 1.0e8, 2.0, 3.0*r_s, 1.0e4*r_s);

        EnergyGrid g(1.0e5, 1.0e11, 60);
        ConstantY xs(1.0e-33, 5.0e-34, 1.0e6, 0.36, 0.36);
        CascadeKernel k(g, xs, 1.0e-12);

        std::vector<double> phi0(g.n());
        for (int i=0;i<g.n();++i) phi0[i] = g.powerlaw_bin_integral(i, 2.0);

        IntegratorOpts o; o.ode_rel_tol = 1.0e-12;
        std::printf("    %-10s %-13s %-11s %-13s %s\n",
                    "b/r_s","E_loc/E_inf","tau","dif em tau","dif em P");
        double pior = 0.0;
        for (double frac : {3.0, 10.0, 100.0}) {
            Ray ray; ray.E_inf_GeV = 1.0e8; ray.b_cm = frac*r_s;
            const double rt = schw.r_turning(ray.b_cm, ray.b_cm);
            const double boost = 1.0/std::sqrt(schw.f(rt));

            const CascadeResult a = transport_cascade(phi0,k,ray,schw,halo,o,false);
            const CascadeResult c = transport_cascade(phi0,k,ray,schw,halo,o,true);

            double sa=0, sc=0, s0=0;
            for (int i=0;i<g.n();++i){ sa+=a.phi[i]; sc+=c.phi[i]; s0+=phi0[i]; }

            // O observavel certo e a diferenca em TAU (equivalente a
            // sigma_eff), nao em P. O limite sqrt(3) controla sigma;
            // P = exp(-tau) amplifica qualquer mudanca de sigma por um
            // fator exp(delta*tau), entao medir P mistura o efeito
            // fisico com a espessura optica escolhida.
            const double tau_a = -std::log(sa/s0);
            const double tau_c = -std::log(sc/s0);
            const double d_tau = std::fabs(tau_a - tau_c)/tau_c;
            const double d_P   = std::fabs(sa - sc)/sc;
            pior = std::max(pior, d_tau);
            std::printf("    %-10.1f %-13.5f %-11.3f %-13.4f %.4f\n",
                        frac, boost, tau_c, d_tau, d_P);
        }
        // Limite: E_loc/E_inf <= sqrt(3), logo com sigma ~ E^0.36 a razao
        // de secoes de choque e no maximo 3^0.18 = 1.216. Esse e o teto
        // ABSOLUTO da diferenca em tau, atingido so por um raio que
        // ficasse na esfera de fotons o caminho inteiro.
        const double teto = std::pow(3.0, 0.18) - 1.0;
        std::printf("    limite: E_loc/E_inf <= sqrt(3), logo dif em tau <= 3^0.18 - 1 = %.4f\n", teto);
        ++total;
        if (!(pior <= teto)) { ++falhas;
            std::printf("    <-- FALHA: dif em tau acima do teto fisico\n"); }
        else std::printf("    OK: dif maxima em tau = %.2f%%, dentro do teto de %.2f%%\n",
                         100*pior, 100*teto);

        // Corolario que vale registrar: a diferenca em P NAO e limitada
        // por isso. P = exp(-tau) amplifica por exp(delta*tau), entao em
        // tau ~ 10 uma diferenca de 4% em sigma vira 32% em P. Na janela
        // observavel (tau em [0.3, 5]) fica em poucos por cento, mas a
        // afirmacao "o redshift nao importa" so vale ali.
        std::printf("    nota: a diferenca em P nao herda esse teto -- ela e\n");
        std::printf("          amplificada por exp(delta*tau). Ver a coluna acima.\n");
    }

    // =================================================================
    std::printf("\n  ROUND-TRIP  do leitor de tabela, sem depender do dipole\n");
    {
        // Gera tabela a partir de kernel ANALITICO, le de volta, e exige
        // reproduzir o kernel. Valida o LEITOR sozinho -- a tabela real
        // de dsigma_NC/dy do dipole ainda nao existe.
        const std::string arq = "/tmp/phasis_rt_dsigma.dat";
        PowerLawY orig(1.0e-33, 5.0e-34, 1.0e6, 0.36, 0.0, 1.0, 1.0e-3);

        write_dsigma_table(arq, orig, 1.0e5, 1.0e11, 25, 1.0e-3, 25, "teste round-trip");

        auto cc = std::make_shared<PowerLawCrossSection>(5.0e-34, 1.0e6, 0.0);
        TableDifferentialCrossSection::Expect esp;
        TableDifferentialCrossSection lida(arq, esp, cc);

        // Uma lei de potencia e reta em log-log, entao bilinear em
        // log-log a reproduz EXATAMENTE, inclusive fora dos nos.
        double pior = 0.0; int n = 0;
        for (int a = 0; a < 17; ++a) {
            const double E = 1.0e5*std::pow(1.0e6, (a + 0.5)/17.0);
            for (int b = 0; b < 19; ++b) {
                const double y = 1.0e-3*std::pow(1.0e3, (b + 0.5)/19.0);
                const double x1 = orig.dsigma_nc_dy(E, y);
                const double x2 = lida.dsigma_nc_dy(E, y);
                if (x1 > 0.0) { pior = std::max(pior, std::fabs(x2-x1)/x1); ++n; }
            }
        }
        std::printf("    %d pontos fora dos nos, pior desvio %.3e\n", n, pior);
        chk("round-trip reproduz o kernel analitico", pior, 0.0, 1.0e-10);

        // Metadados: ausencia de QUALQUER chave tem de lancar.
        int lancou = 0;
        const char* chaves[] = {"convention_y","target","projectile","current",
                                "units_sigma","units_E","M_Z_GeV","dipole_model","generated_by"};
        for (const char* ch : chaves) {
            const std::string mut = "/tmp/phasis_rt_mut.dat";
            std::ifstream in(arq); std::ofstream out(mut);
            std::string ln;
            while (std::getline(in, ln)) {
                if (ln.size() > 1 && ln[0] == '#' && ln.find(ch) != std::string::npos) continue;
                out << ln << "\n";
            }
            out.close();
            try { TableDifferentialCrossSection t(mut, esp, cc); }
            catch (const std::runtime_error&) { ++lancou; }
        }
        std::printf("    %d de 9 chaves obrigatorias fazem o leitor lancar quando ausentes\n", lancou);
        chk("todo metadado e obrigatorio", static_cast<double>(lancou), 9.0, 0.0);

        // Convencao divergente tem de lancar, nao ser aceita em silencio.
        bool rec = false;
        TableDifferentialCrossSection::Expect outra;
        outra.convention_y = "(E_in - E_out)/E_out";
        try { TableDifferentialCrossSection t(arq, outra, cc); }
        catch (const std::runtime_error&) { rec = true; }
        chk("convencao divergente lanca", rec ? 1.0 : 0.0, 1.0, 0.0);
    }

    std::printf("\n=====================================================\n");
    std::printf("  %d verificacoes, %d falhas\n", total, falhas);
    std::printf("  RESULTADO: %s\n", falhas==0?"OK":"FALHA");
    std::printf("=====================================================\n\n");
    return falhas==0?0:1;
}
