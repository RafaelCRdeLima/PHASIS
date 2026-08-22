#ifndef PHASIS_SWEEP_HPP
#define PHASIS_SWEEP_HPP

#include <string>
#include <vector>

#include "phasis/trace.hpp"

namespace phasis {

// =====================================================================
// Varredura sobre muitas geodesicas.
//
// trace_ray ja e funcao pura, entao o paralelismo aqui e um
// `#pragma omp parallel for` sobre um vector<Ray> plano, escrevendo em
// posicoes pre-alocadas de vector<Result>. Sem RNG, sem estado mutavel
// compartilhado, sem I/O dentro do laco -- o que torna o resultado
// determinista e independente do numero de threads (T18).
//
// Duas coisas que NAO podem faltar:
//
//  1. EXCECOES. TableCrossSection lanca fora da faixa tabelada. Uma
//     excecao atravessando a fronteira de uma regiao OpenMP e
//     comportamento indefinido. O corpo do laco fica dentro de
//     try/catch e converte a falha em Result{error, error_msg}. Nada
//     escapa.
//
//  2. AGREGACAO. max_evals e tolerance_met sao por-raio. Com 1e5 raios
//     ninguem le linha por linha, entao o sumario conta quantos
//     degradaram, quantos sao quase-criticos, quantos foram capturados,
//     e guarda uma amostra das mensagens de erro.
// =====================================================================

struct SweepSummary {
    long n_rays          = 0;
    long n_captured      = 0;
    long n_near_critical = 0;
    long n_degraded      = 0;   // tolerance_met == false
    long n_errors        = 0;
    long n_no_matter     = 0;

    double tau_min  = 0.0;
    double tau_max  = 0.0;
    double tau_mean = 0.0;      // sobre os raios finitos e sem erro

    long n_evals_total = 0;

    // Amostra das mensagens distintas, para nao inundar a saida.
    std::vector<std::string> error_sample;

    std::string to_string() const;
};

// results e redimensionado para rays.size() e preenchido em paralelo.
SweepSummary sweep(const std::vector<Ray>& rays,
                   const Metric& metric,
                   const DensityProfile& profile,
                   const CrossSection& xsec,
                   const IntegratorOpts& opts,
                   std::vector<Result>& results);

} // namespace phasis

#endif
