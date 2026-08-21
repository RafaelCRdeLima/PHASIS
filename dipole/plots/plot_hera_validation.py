#!/usr/bin/env python3
"""F1 - sigma_red do HERA: dados vs limite EM do modelo de dipolo.

    python3 plots/plot_hera_validation.py [arquivo.dat ...]

Sem argumento usa data/hera_validation.dat. Com varios, sobrepoe as curvas
(util para comparar grades de quadratura).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def carrega(path):
    d = np.loadtxt(path)
    meta = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith("#"):
            break
        if "=" in line:
            k, _, v = line[1:].partition("=")
            meta[k.strip()] = v.strip()
    return d, meta


def main(argv):
    arquivos = argv[1:] or ["data/hera_validation.dat"]
    dados, metas = zip(*(carrega(a) for a in arquivos))

    # bins de Q2 para os paineis
    Q2_todos = np.unique(np.round(dados[0][:, 0], 6))
    escolhidos = Q2_todos[np.linspace(0, len(Q2_todos) - 1, min(6, len(Q2_todos))).astype(int)]

    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), sharey=True)
    cores = plt.cm.viridis(np.linspace(0.15, 0.8, len(arquivos)))

    for ax, Q2 in zip(axes.ravel(), escolhidos):
        d0 = dados[0]
        m = np.isclose(d0[:, 0], Q2)
        ax.errorbar(d0[m, 1], d0[m, 3], yerr=d0[m, 4], fmt="o", ms=4,
                    color="k", label="HERA I+II", zorder=3)
        for d, meta, cor, nome in zip(dados, metas, cores, arquivos):
            mm = np.isclose(d[:, 0], Q2)
            o = np.argsort(d[mm, 1])
            rot = meta.get("quadratura", Path(nome).stem)
            ax.plot(d[mm, 1][o], d[mm, 5][o], "-", lw=1.8, color=cor, label=rot)
        ax.set_xscale("log")
        ax.set_title(r"$Q^2 = %g$ GeV$^2$" % Q2, fontsize=10)
        ax.set_xlabel("$x$")
        ax.grid(alpha=0.25, lw=0.5)

    axes[0, 0].set_ylabel(r"$\sigma_{\rm red}$")
    axes[1, 0].set_ylabel(r"$\sigma_{\rm red}$")
    axes[0, 0].legend(fontsize=7.5, loc="best")
    fig.suptitle(r"F1 — limite eletromagnético vs $\sigma_{\rm red}$ do HERA "
                 r"(H1+ZEUS, EPJ C75 (2015) 580)", fontsize=12)
    fig.tight_layout()
    out = "plots/hera_validation.png"
    fig.savefig(out, dpi=140)
    print("Figura gerada:", out)

    # painel de pull
    fig2, ax = plt.subplots(figsize=(9, 4.2))
    for d, meta, cor, nome in zip(dados, metas, cores, arquivos):
        rot = meta.get("quadratura", Path(nome).stem)
        chi2 = np.mean(d[:, 8] ** 2)
        ax.plot(d[:, 0], d[:, 8], "o", ms=3.5, color=cor, alpha=0.7,
                label=r"%s   $\chi^2$/pt = %.4g" % (rot, chi2))
    ax.axhline(0, color="k", lw=0.8)
    for s in (-1, 1):
        ax.axhline(s, color="k", lw=0.6, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel(r"$Q^2$ [GeV$^2$]")
    ax.set_ylabel("pull  (modelo − dado)/erro")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25, lw=0.5)
    fig2.tight_layout()
    out2 = "plots/hera_validation_pulls.png"
    fig2.savefig(out2, dpi=140)
    print("Figura gerada:", out2)


if __name__ == "__main__":
    main(sys.argv)
