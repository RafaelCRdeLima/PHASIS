#include "phasis/sweep.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <sstream>

namespace phasis {

SweepSummary sweep(const std::vector<Ray>& rays,
                   const Metric& metric,
                   const DensityProfile& profile,
                   const CrossSection& xsec,
                   const IntegratorOpts& opts,
                   std::vector<Result>& results)
{
    const long n = static_cast<long>(rays.size());
    results.assign(rays.size(), Result{});

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
    for (long i = 0; i < n; ++i) {
        try {
            results[static_cast<std::size_t>(i)] =
                trace_ray(rays[static_cast<std::size_t>(i)],
                          metric, profile, xsec, opts);
        }
        catch (const std::exception& e) {
            // Nada pode escapar da regiao paralela: excecao atravessando
            // a fronteira de OpenMP e comportamento indefinido.
            Result r;
            r.error = true;
            r.error_msg = e.what();
            results[static_cast<std::size_t>(i)] = r;
        }
        catch (...) {
            Result r;
            r.error = true;
            r.error_msg = "excecao desconhecida";
            results[static_cast<std::size_t>(i)] = r;
        }
    }

    // Agregacao em serie, depois do laco: ordem fixa, resultado
    // independente do numero de threads.
    SweepSummary s;
    s.n_rays = n;

    bool primeiro = true;
    double soma = 0.0;
    long n_finitos = 0;

    for (const Result& r : results) {
        s.n_evals_total += r.n_evals;

        if (r.error) {
            ++s.n_errors;
            if (std::find(s.error_sample.begin(), s.error_sample.end(), r.error_msg)
                == s.error_sample.end() && s.error_sample.size() < 5) {
                s.error_sample.push_back(r.error_msg);
            }
            continue;
        }
        if (r.captured)         ++s.n_captured;
        if (r.near_critical)    ++s.n_near_critical;
        if (!r.tolerance_met)   ++s.n_degraded;
        if (!r.crosses_matter)  ++s.n_no_matter;

        if (std::isfinite(r.tau)) {
            soma += r.tau;
            ++n_finitos;
            if (primeiro) { s.tau_min = s.tau_max = r.tau; primeiro = false; }
            else {
                s.tau_min = std::min(s.tau_min, r.tau);
                s.tau_max = std::max(s.tau_max, r.tau);
            }
        }
    }
    if (n_finitos > 0) s.tau_mean = soma/static_cast<double>(n_finitos);

    return s;
}

std::string SweepSummary::to_string() const
{
    std::ostringstream m;
    m << "raios            : " << n_rays << "\n"
      << "  capturados     : " << n_captured << "\n"
      << "  quase-criticos : " << n_near_critical << "\n"
      << "  sem materia    : " << n_no_matter << "\n"
      << "  degradados     : " << n_degraded
      << "   (tolerancia nao atingida)\n"
      << "  com erro       : " << n_errors << "\n";
    for (const std::string& e : error_sample) m << "      * " << e << "\n";
    m << "tau (finitos)    : min " << tau_min
      << "  media " << tau_mean << "  max " << tau_max << "\n"
      << "avaliacoes       : " << n_evals_total << "\n";
    return m.str();
}

} // namespace phasis
