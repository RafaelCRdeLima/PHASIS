#ifndef PHASIS_ODE_HPP
#define PHASIS_ODE_HPP

#include <array>
#include <cstddef>
#include <functional>
#include <vector>

namespace phasis {

// =====================================================================
// Dormand-Prince 5(4) com passo adaptativo.
//
// Existe porque, com densidade dependente de theta, psi(r) e necessario
// DENTRO do integrando de tau: as quadraturas deixam de ser
// independentes. Resolver isso com sub-quadratura aninhada custaria
// O(N^2) por raio, inviavel para 1e5 raios. Aqui o estado (psi, l, tau)
// avanca junto, numa unica passada.
//
// Sem dense output: o acoplamento e resolvido levando tudo no mesmo
// vetor de estado, entao nunca e preciso interpolar psi a posteriori.
// =====================================================================

template <std::size_t N>
using OdeState = std::array<double, N>;

struct OdeOpts {
    double rel_tol   = 1.0e-12;

    // Piso ABSOLUTO do controle de erro. Nao pode ser ~0.
    //
    // O controlador usa sc_i = abs_tol + rel_tol*|y_i|. Se abs_tol for
    // desprezivel e alguma componente valer, digamos, 1e-80 (o ramo do
    // raio que se afasta do disco tem tau ~ 1e-80), entao
    // sc_i ~ 1e-89 e qualquer ruido nessa componente da erro gigante: o
    // passo encolhe ate h_min por causa de uma quantidade fisicamente
    // irrelevante, e o resultado sai marcado como degradado sem que
    // nada de errado tenha acontecido.
    //
    // 1e-10 e desprezivel para TODAS as componentes do estado deste
    // problema (tau ~ 1e2, psi ~ 1 rad, l ~ 1e9 cm, X ~ 1e11 g/cm^2) e
    // elimina o colapso. Foi diagnosticado assim: 180 de 10000 raios em
    // T18 saiam degradados, todos com b logo acima de b_crit, com poucos
    // passos -- ou seja h < h_min, nao max_steps.
    double abs_tol   = 1.0e-10;
    double h_init    = 0.0;     // 0 = escolher automaticamente
    double h_min_rel = 1.0e-14; // relativo ao intervalo total

    // Teto de passo, relativo ao intervalo total. Nao e cosmetico: num
    // trecho onde o integrando e identicamente nulo o controlador nao ve
    // erro algum e cresce o passo 5x por iteracao, podendo pular por
    // cima de estrutura fina (a espessura vertical de um disco, por
    // exemplo). O teto limita o estrago; a defesa principal e cortar o
    // intervalo nas fronteiras do suporte, o que trace_via_ode faz.
    double h_max_rel = 0.05;

    // Piso absoluto POR COMPONENTE. Se vazio, usa o escalar abs_tol.
    //
    // Necessario na cascata: os bins de alta energia caem muitas ordens
    // enquanto os de baixa crescem, no MESMO vetor de estado. Um piso
    // global escalado pelo maior componente deixa os pequenos sem
    // controle; escalado pelo menor, colapsa o passo (o bug da Fase 3).
    // A escolha certa e um piso por bin, proporcional ao fluxo INICIAL
    // daquele bin: diz "nao me importo com precisao relativa depois que
    // este bin caiu muito abaixo do que ele proprio comecou", que e
    // fisicamente o certo -- um bin nessa situacao ja nao contribui.
    std::vector<double> abs_tol_per_component;
    long   max_steps = 2000000;
};

struct OdeStats {
    long steps    = 0;
    long rejected = 0;
    long n_evals  = 0;
    bool ok       = true;       // false se hit_min_step ou hit_max_steps

    // Distinguir os dois modos importa: h_min significa que o
    // controlador nao conseguiu satisfazer a tolerancia (estrutura fina
    // demais, ou tolerancia pedida abaixo do que a aritmetica permite),
    // enquanto max_steps significa so que o intervalo e longo demais
    // para o orcamento. As acoes corretivas sao opostas.
    bool hit_min_step  = false;
    bool hit_max_steps = false;
    double h_final     = 0.0;
};

// Integra dy/ds = f(s, y) de s0 a s1. y entra com a condicao inicial e
// sai com o resultado.
template <std::size_t N>
OdeStats dopri54(const std::function<void(double, const OdeState<N>&, OdeState<N>&)>& f,
                 double s0, double s1, OdeState<N>& y, const OdeOpts& opts);

// instanciacao usada pelo tracador
// Versao de dimensao dinamica, para a cascata (N bins nao e conhecido em
// tempo de compilacao).
OdeStats dopri54_dyn(
    const std::function<void(double, const std::vector<double>&, std::vector<double>&)>& f,
    double s0, double s1, std::vector<double>& y, const OdeOpts& opts);

extern template OdeStats dopri54<6>(
    const std::function<void(double, const OdeState<6>&, OdeState<6>&)>&,
    double, double, OdeState<6>&, const OdeOpts&);

} // namespace phasis

#endif
