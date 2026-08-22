#ifndef PHASIS_SHELL_HPP
#define PHASIS_SHELL_HPP

#include <functional>
#include <string>
#include <vector>

namespace phasis {

// =====================================================================
// Localizacao da casca tau = alvo, ao longo de um parametro escalar.
//
// Por que existe: com ~25 ordens de grandeza de variacao em tau, uma
// grade uniforme gasta quase tudo onde a resposta e 0 ou 1. A fisica
// vive na casca fina onde tau ~ 1, e e la que a sensibilidade ao modelo
// de sigma existe -- fora de tau em [0.3, 5] o resultado e insensivel
// ao modelo, por melhor que ele seja.
//
// ARMADILHA PRINCIPAL, e a razao do protocolo abaixo: bissecao assume
// monotonicidade, e tau(parametro) NAO e monotonico em geral. Um perfil
// oco ja basta para quebrar: a corda por uma casca [r1,r2] CRESCE com b
// ate b = r1 (de 2(r2-r1) para 2*sqrt(r2^2-r1^2)) e so depois cai. Com
// disco flared e raios encurvados e pior. Bissecao cega devolveria UMA
// raiz em silencio.
//
// Protocolo, nesta ordem:
//   1. varredura grosseira (>= 32 pontos) contando mudancas de sinal;
//   2. exatamente uma  -> bissecta dentro do bracket;
//   3. mais de uma     -> registra TODAS e marca multi_root;
//   4. nenhuma         -> shell_found = false com o valor extremo.
//                         Sem extrapolar, sem iterar ate max_iter, sem NaN.
// =====================================================================

struct ShellOpts {
    double target_tau = 1.0;
    int    n_coarse   = 64;     // minimo util: 32
    double rel_tol    = 1.0e-12;
    int    max_iter   = 200;
};

struct ShellResult {
    bool shell_found = false;
    bool multi_root  = false;

    std::vector<double> roots;      // todos os valores do parametro com tau = alvo

    double tau_min = 0.0;           // extremos vistos na varredura grosseira
    double tau_max = 0.0;
    double p_at_tau_min = 0.0;
    double p_at_tau_max = 0.0;

    int  n_coarse       = 0;
    int  n_sign_changes = 0;
    long n_tau_evals    = 0;

    std::string diagnostico() const;
};

// tau_of(p) deve devolver tau >= 0. Valores 0 e +inf sao aceitos e tem
// sinal bem definido em log; NaN e erro do chamador.
ShellResult locate_shell(const std::function<double(double)>& tau_of,
                         double p_lo, double p_hi,
                         const ShellOpts& opts = ShellOpts{});

} // namespace phasis

#endif
