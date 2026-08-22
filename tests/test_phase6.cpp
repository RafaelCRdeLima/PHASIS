// =====================================================================
// Refatoracao da topologia -- T31, T32 (geometria pura).
//
// T31 e o unico teste do projeto cujo valor de referencia e exato,
// geometrico, e independente de tudo o que foi construido: se classify()
// erra, erra por 45 graus, nao por 1e-13.
// =====================================================================

#include "phasis/emission.hpp"
#include "phasis/cross_section.hpp"
#include "phasis/density.hpp"
#include "phasis/metric.hpp"
#include "phasis/trace.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace phasis;

namespace {

constexpr double kPi = 3.141592653589793;
int falhas = 0, total = 0;

void chk(const char* n, double g, double r, double t)
{
    ++total;
    const double d = (std::fabs(r) > 0.0) ? std::fabs(r) : 1.0;
    const double rel = std::fabs(g - r)/d;
    const bool ok = (rel <= t) && std::isfinite(g);
    if (!ok) ++falhas;
    std::printf("    %-46s %-18.12g ref %-18.12g rel %8.2e  %s\n",
                n, g, r, rel, ok ? "OK" : "<-- FALHA");
}

void chk_bool(const char* n, bool g, bool r)
{
    ++total;
    const bool ok = (g == r);
    if (!ok) ++falhas;
    std::printf("    %-46s %-18s ref %-18s                %s\n",
                n, g?"true":"false", r?"true":"false", ok?"OK":"<-- FALHA");
}

} // namespace

int main()
{
    const double r_s = 1.0e6;
    const Schwarzschild schw(r_s);
    const Minkowski flat;
    IntegratorOpts o; o.rel_tol = 1.0e-13;

    std::printf("\n=====================================================\n");
    std::printf(" PHASIS -- topologia do raio emitido (T31, T32)\n");
    std::printf("=====================================================\n");

    // =================================================================
    std::printf("\n  T31a  fronteira de classify() contra sin(psi_c) = b_crit sqrt(f)/r\n");
    {
        std::printf("    %-10s %-16s %-16s %s\n","r/r_s","psi_c (formula)","fronteira varrida","dif rel");
        double pior = 0.0;
        for (double fr : {1.5, 2.0, 3.0, 5.0, 10.0, 100.0}) {
            const double re = fr*r_s;
            const double psi_c = escape_cone_angle(re, schw);

            // A fronteira de captura em psi esta em pi - psi_c (para
            // r > r_ph). Bissecta a mudanca de classificacao.
            double lo = kPi/2, hi = kPi;
            for (int i = 0; i < 200; ++i) {
                const double m = 0.5*(lo + hi);
                if (classify(re, m, schw) == Topology::Captured) hi = m; else lo = m;
            }
            const double fronteira = 0.5*(lo + hi);
            const double esperado  = kPi - psi_c;
            const double d = std::fabs(fronteira - esperado)/esperado;
            pior = std::max(pior, d);
            std::printf("    %-10.1f %-16.12f %-16.12f %.2e\n", fr, psi_c, fronteira, d);
        }
        chk("fronteira coincide com pi - psi_c", pior, 0.0, 1.0e-12);

        // Os dois valores EXATOS
        chk("psi_c(1.5 r_s) = pi/2",  escape_cone_angle(1.5*r_s, schw), kPi/2, 1.0e-14);
        chk("psi_c(3.0 r_s) = pi/4",  escape_cone_angle(3.0*r_s, schw), kPi/4, 1.0e-12);
        chk("f_cap(1.5 r_s) = 1/2",   captured_fraction(1.5*r_s, schw), 0.5, 1.0e-14);
        chk("f_cap(3.0 r_s) = (2-sqrt2)/4",
            captured_fraction(3.0*r_s, schw), (2.0-std::sqrt(2.0))/4.0, 1.0e-12);
    }

    // =================================================================
    std::printf("\n  T31b  fracao capturada por Monte Carlo, cos(psi) uniforme\n");
    {
        // Emissao isotropica no referencial local: cos(psi) uniforme em
        // [-1,1]. Amostrar psi uniforme seria errado.
        const long N = 1000000;
        std::printf("    %-10s %-16s %-16s %-12s %s\n",
                    "r/r_s","f_cap exata","f_cap Monte Carlo","desvio","em sigmas");
        double pior_sig = 0.0;
        for (double fr : {1.5, 2.0, 3.0, 5.0, 10.0, 100.0}) {
            const double re = fr*r_s;
            const double exata = captured_fraction(re, schw);

            long cap = 0;
            for (long k = 0; k < N; ++k) {
                // grade determinista em cos(psi): sem RNG, reprodutivel
                const double c = -1.0 + 2.0*(k + 0.5)/N;
                const double psi = std::acos(c);
                if (classify(re, psi, schw) == Topology::Captured) ++cap;
            }
            const double mc = static_cast<double>(cap)/N;
            const double sig = std::sqrt(std::max(exata*(1.0-exata), 1.0e-12)/N);
            const double nsig = std::fabs(mc - exata)/sig;
            pior_sig = std::max(pior_sig, nsig);
            std::printf("    %-10.1f %-16.12f %-16.12f %-12.2e %.2f\n",
                        fr, exata, mc, std::fabs(mc-exata), nsig);
        }
        // grade determinista: o erro e de discretizacao, nao estatistico,
        // e vale ~1/N. Cobramos 3 sigma binomiais como pedido.
        ++total;
        if (!(pior_sig < 3.0)) { ++falhas;
            std::printf("    <-- FALHA: desvio de %.2f sigma\n", pior_sig); }
        else std::printf("    OK: pior desvio %.2f sigma (limite 3)\n", pior_sig);
    }

    // =================================================================
    std::printf("\n  T31c  a QUARTA linha: abaixo da esfera de fotons a condicao INVERTE\n");
    {
        const double re = 1.2*r_s;                 // dentro da esfera de fotons
        const double b_crit = schw.b_crit();
        std::printf("    r_emit = 1.2 r_s < 1.5 r_s ;  b_crit = %.6e\n", b_crit);
        std::printf("    %-12s %-14s %-14s %s\n","psi","b/b_crit","topologia","esperado");

        int erros = 0;
        for (double psi : {0.05, 0.3, 0.7, 1.0, 1.5, 2.0, 3.0}) {
            double b = 0.0;
            const Topology t = classify(re, psi, schw, &b);
            const bool fora = (psi < kPi/2);
            const Topology esp = fora
                ? ((b < b_crit) ? Topology::OutboundOnly : Topology::Captured)
                : Topology::Captured;
            if (t != esp) ++erros;
            std::printf("    %-12.2f %-14.6f %-14s %s\n",
                        psi, b/b_crit, to_string(t), to_string(esp));
        }
        ++total; if (erros) ++falhas;
        std::printf("    %s\n", erros ? "<-- FALHA" : "OK: condicao invertida reproduzida");

        // Contraste com r_emit > r_ph, onde a condicao e a normal
        const double re2 = 3.0*r_s;
        double b2 = 0.0;
        chk_bool("r>r_ph, psi=0.3 (fora): OutboundOnly",
                 classify(re2, 0.3, schw, &b2) == Topology::OutboundOnly, true);
        chk_bool("r<r_ph, psi=0.3 (fora): depende de b",
                 classify(re, 0.3, schw, &b2) != Topology::Turning, true);
    }

    // =================================================================
    std::printf("\n  T32  continuidade na esfera de fotons\n");
    {
        std::printf("    %-16s %-18s %-18s %s\n",
                    "r_emit/r_s","f_cap","psi_c","classify(psi=pi/2 -+ eps)");
        for (double d : {1.0e-3, 1.0e-6, 1.0e-9}) {
            for (int lado = 0; lado < 2; ++lado) {
                const double re = 1.5*r_s*(lado ? (1.0 + d) : (1.0 - d));
                const double fc = captured_fraction(re, schw);
                const double pc = escape_cone_angle(re, schw);
                double bb = 0.0;
                const Topology tm = classify(re, kPi/2 - 1.0e-12, schw, &bb);
                const Topology tp = classify(re, kPi/2 + 1.0e-12, schw, &bb);
                std::printf("    %-16.12f %-18.12f %-18.12f %s / %s\n",
                            re/r_s, fc, pc, to_string(tm), to_string(tp));
                ++total;
                if (std::fabs(fc - 0.5) > 1.0e-3) { ++falhas;
                    std::printf("      <-- FALHA: f_cap longe de 1/2\n"); }
            }
        }
        std::printf("    (f_cap -> 1/2 pelos DOIS lados: a classificacao e continua\n");
        std::printf("     na fronteira, e o exemplo_buraco_negro.cfg poe r_in exatamente ali)\n");

        // r_in do config de estresse, exatamente na fronteira
        chk("f_cap em r = 1.5 r_s exato", captured_fraction(1.5*r_s, schw), 0.5, 1.0e-14);
    }

    std::printf("\n  T33  consistencia com a Fase 2 (infinito a infinito)\n");
    {
        // Um raio Turning emitido de r_emit MUITO grande, para dentro,
        // reproduz o resultado infinito-a-infinito para o mesmo b,
        // descontada a coluna alem de r_emit. Com r_emit >= r_out nao ha
        // materia alem, entao os dois tem de coincidir.
        const PowerLawHalo halo(1.0e2, 1.0e8, 2.0, 5.0*r_s, 2.0e3*r_s);
        const PowerLawCrossSection xs(1.0e-33, 1.0e6, 0.36);

        std::printf("    %-10s %-18s %-18s %s\n","b/r_s","Fase 2 (inf-inf)","emitido r_emit=r_out","dif rel");
        double pior=0;
        for (double fr : {4.0, 20.0, 200.0}) {
            Ray a; a.E_inf_GeV=1.0e9; a.b_cm=fr*r_s;              // do infinito
            const double ta = trace_ray(a,schw,halo,xs,o).tau;

            Ray b = a;                                            // emitido na borda
            b.r_emit_cm = halo.r_support_max();
            b.outward   = false;                                  // para dentro
            const double tb = trace_ray(b,schw,halo,xs,o).tau;

            const double d=std::fabs(tb-ta)/ta; pior=std::max(pior,d);
            std::printf("    %-10.1f %-18.12e %-18.12e %.2e\n",fr,ta,tb,d);
        }
        chk("emitido em r_out == vindo do infinito",pior,0.0,1.0e-11);

        // Emitido no MEIO do halo, para dentro: tem de dar menos que o
        // do infinito, e a diferenca e exatamente o ramo de entrada
        // acima de r_emit.
        Ray c; c.E_inf_GeV=1.0e9; c.b_cm=20.0*r_s;
        const Result inf = trace_ray(c,schw,halo,xs,o);
        Ray d2=c; d2.r_emit_cm=3.0e2*r_s; d2.outward=false;
        const Result em = trace_ray(d2,schw,halo,xs,o);
        std::printf("    emitido em r=300 r_s (dentro): tau=%.6e vs %.6e do infinito\n",
                    em.tau, inf.tau);
        chk("saida identica; entrada truncada",em.tau_outbound,inf.tau_outbound,1.0e-11);
        ++total; if(!(em.tau_inbound < inf.tau_inbound)){++falhas;
            std::printf("    <-- FALHA: entrada nao foi truncada\n");}
        else std::printf("    OK: entrada %.4e < %.4e\n",em.tau_inbound,inf.tau_inbound);
    }

    std::printf("\n  T34  radial exato (b = 0) vs caminho geral com b pequeno\n");
    {
        const UniformBall bola(2.0, 1.0e9);
        const PowerLawCrossSection xs(1.0e-33, 1.0e6, 0.36);

        // Minkowski, do infinito: nada captura.
        {
            Ray a; a.E_inf_GeV=1.0e9; a.b_cm=0.0;
            Ray b=a; b.b_cm=1.0e-8*bola.r_support_max();
            chk("Minkowski, do infinito: b=0 vs b=1e-8 R",
                trace_ray(a,flat,bola,xs,o).tau, trace_ray(b,flat,bola,xs,o).tau, 1.0e-10);
        }
        // Schwarzschild: um raio radial VINDO DO INFINITO cai no buraco
        // negro -- e isso e fisica, nao falha. O caso radial que escapa e
        // o EMITIDO para fora, que so existe com a topologia nova.
        {
            Ray a; a.E_inf_GeV=1.0e9; a.b_cm=0.0;
            chk("Schwarzschild, radial do infinito: capturado",
                trace_ray(a,schw,bola,xs,o).captured ? 1.0 : 0.0, 1.0, 0.0);

            const double re = 20.0*r_s;
            const Ray p0 = emit_ray(re, 0.0,     1.0e9, schw);   // radial para fora
            const Ray p1 = emit_ray(re, 1.0e-8,  1.0e9, schw);   // quase radial
            const Result r0 = trace_ray(p0,schw,bola,xs,o);
            const Result r1 = trace_ray(p1,schw,bola,xs,o);
            std::printf("    emitido para fora de r=20 r_s: b=%.3e e %.3e\n", p0.b_cm, p1.b_cm);
            chk("Schwarzschild, emitido radial: b=0 vs b pequeno", r0.tau, r1.tau, 1.0e-10);
            chk("radial nao capturado", r0.captured?1.0:0.0, 0.0, 0.0);
        }
        // psi = 0 e psi = pi pela via de emissao dao b = 0 exato
        const double be0 = b_from_emission(1.0e8, 0.0, schw);
        const double bep = b_from_emission(1.0e8, 3.141592653589793, schw);
        chk("b(psi=0)",be0,0.0,0.0);
        // sin(M_PI) vale 1.22e-16, nao 0, porque M_PI nao e exatamente
        // pi. Normalizado por r_emit isso e precisao de maquina.
        chk("b(psi=pi)/r_emit",bep/1.0e8,0.0,1.0e-15);
    }

    std::printf("\n=====================================================\n");
    std::printf("  %d verificacoes, %d falhas\n", total, falhas);
    std::printf("  RESULTADO: %s\n", falhas == 0 ? "OK" : "FALHA");
    std::printf("=====================================================\n\n");
    return falhas == 0 ? 0 : 1;
}
