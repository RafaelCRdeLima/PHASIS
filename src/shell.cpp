#include "phasis/shell.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace phasis {

namespace {

// g(p) = ln(tau(p)) - ln(alvo). A raiz e g = 0.
//
// Trabalhar em log nao e cosmetico: tau varre ~25 ordens de grandeza, e
// em escala linear a bissecao gastaria toda a precisao no lado opaco.
//
// tau = 0 -> g = -inf ; tau = +inf (capturado) -> g = +inf. Os dois tem
// SINAL bem definido, que e tudo o que o bracketing precisa.
double g_de(double tau, double alvo)
{
    if (!(tau >= 0.0)) {
        throw std::runtime_error("locate_shell: tau negativo ou NaN");
    }
    if (tau == 0.0)              return -std::numeric_limits<double>::infinity();
    if (std::isinf(tau))         return  std::numeric_limits<double>::infinity();
    return std::log(tau) - std::log(alvo);
}

int sinal(double g)
{
    if (g > 0.0) return  1;
    if (g < 0.0) return -1;
    return 0;
}

} // namespace

ShellResult locate_shell(const std::function<double(double)>& tau_of,
                         double p_lo, double p_hi,
                         const ShellOpts& opts)
{
    ShellResult out;

    if (!(p_hi > p_lo)) throw std::runtime_error("locate_shell: intervalo invalido");
    const int n = std::max(32, opts.n_coarse);
    out.n_coarse = n;

    // ---- 1. varredura grosseira -------------------------------------
    std::vector<double> p(static_cast<std::size_t>(n)), tau(static_cast<std::size_t>(n)),
                        g(static_cast<std::size_t>(n));

    for (int k = 0; k < n; ++k) {
        p[static_cast<std::size_t>(k)] = p_lo + (p_hi - p_lo)*k/(n - 1);
        tau[static_cast<std::size_t>(k)] = tau_of(p[static_cast<std::size_t>(k)]);
        g[static_cast<std::size_t>(k)]   = g_de(tau[static_cast<std::size_t>(k)], opts.target_tau);
        ++out.n_tau_evals;
    }

    out.tau_min = out.tau_max = tau[0];
    out.p_at_tau_min = out.p_at_tau_max = p[0];
    for (int k = 1; k < n; ++k) {
        const std::size_t q = static_cast<std::size_t>(k);
        if (tau[q] < out.tau_min) { out.tau_min = tau[q]; out.p_at_tau_min = p[q]; }
        if (tau[q] > out.tau_max) { out.tau_max = tau[q]; out.p_at_tau_max = p[q]; }
    }

    // ---- 2. conta mudancas de sinal ---------------------------------
    std::vector<std::pair<int,int>> brackets;
    for (int k = 0; k + 1 < n; ++k) {
        const int sa = sinal(g[static_cast<std::size_t>(k)]);
        const int sb = sinal(g[static_cast<std::size_t>(k+1)]);
        if (sa == 0) { out.roots.push_back(p[static_cast<std::size_t>(k)]); continue; }
        if (sa*sb < 0) brackets.emplace_back(k, k+1);
    }
    out.n_sign_changes = static_cast<int>(brackets.size());

    // ---- 4. nenhum bracket: nao ha casca ----------------------------
    // Devolve limpo, com os extremos. NAO extrapola, NAO itera ate
    // max_iter, NAO devolve NaN.
    if (brackets.empty() && out.roots.empty()) {
        out.shell_found = false;
        return out;
    }

    // ---- 3. bissecta CADA bracket -----------------------------------
    // Todos, nao so o primeiro: com mais de uma raiz, escolher uma em
    // silencio e o modo de falha que nao produz erro visivel.
    for (const auto& br : brackets) {
        double a = p[static_cast<std::size_t>(br.first)];
        double b = p[static_cast<std::size_t>(br.second)];
        double ga = g[static_cast<std::size_t>(br.first)];

        for (int it = 0; it < opts.max_iter; ++it) {
            const double m  = 0.5*(a + b);
            const double gm = g_de(tau_of(m), opts.target_tau);
            ++out.n_tau_evals;

            if (gm == 0.0) { a = b = m; break; }
            if (sinal(ga)*sinal(gm) < 0) { b = m; }
            else                          { a = m; ga = gm; }

            if (std::fabs(b - a) <= opts.rel_tol*std::max(std::fabs(a), std::fabs(b))) break;
        }
        out.roots.push_back(0.5*(a + b));
    }

    std::sort(out.roots.begin(), out.roots.end());
    out.shell_found = !out.roots.empty();
    out.multi_root  = (out.roots.size() > 1);
    return out;
}

std::string ShellResult::diagnostico() const
{
    std::ostringstream m;
    if (!shell_found) {
        m << "casca NAO encontrada: tau em ["
          << tau_min << ", " << tau_max << "] em toda a faixa varrida"
          << (tau_max < 1.0 ? "  (todo transparente)" : "")
          << (tau_min > 1.0 ? "  (todo opaco)" : "");
        return m.str();
    }
    m << roots.size() << " raiz(es)";
    if (multi_root) m << "  [MULTI_ROOT: tau nao e monotonico neste intervalo]";
    m << " em p =";
    for (double r : roots) m << " " << r;
    m << ";  " << n_sign_changes << " mudancas de sinal em "
      << n_coarse << " pontos grosseiros;  " << n_tau_evals << " avaliacoes";
    return m.str();
}

} // namespace phasis
