"""Comparacao com a literatura.

REFERENCIA UNICA, e ela esta em refs/ neste repositorio:

  Formaggio & Zeller, "From eV to EeV: Neutrino Cross-Sections Across
  Energy Scales", Rev. Mod. Phys. 84 (2012) 1307 [arXiv:1305.7513],
  eqs. (90) e (91), citando Gandhi et al. (1996):

      sigma_CC = 5.53e-36 cm^2 (E/GeV)^alpha
      sigma_NC = 2.31e-36 cm^2 (E/GeV)^alpha       alpha ~= 0.363

  validas para 1e16 eV <= E <= 1e21 eV, isto e 1e7 a 1e12 GeV.

As curvas de GQRS sao desenhadas SO nessa faixa. Estende-las seria
atribuir a referencia uma afirmacao que ela nao faz -- e justamente na
regiao onde a saturacao muda a resposta.

Kutak-Kwiecinski EPJ C29 (2003) 521, que este codigo implementa, esta em
refs/s2003-01236-y.pdf mas NAO tabela sigma: so tem figuras. Comparar
com ela exigiria digitalizar curvas, o que nao se faz aqui.
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

mpl.rcParams.update({
    "font.family": "serif", "font.size": 9.5,
    "axes.linewidth": 0.7, "axes.labelsize": 10,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "legend.frameon": False, "legend.fontsize": 8.2,
    "figure.dpi": 140, "savefig.bbox": "tight",
})

C = {"gbw": "#1f4e9c", "iim": "#c85a1e", "col": "#2e7d4f",
     "lit": "#111111", "aux": "#9aa0a8"}

# --- Gandhi et al. via Formaggio & Zeller eqs. (90), (91) ------------
GQRS_ALPHA = 0.363
GQRS_CC, GQRS_NC = 5.53e-36, 2.31e-36
GQRS_LO, GQRS_HI = 1.0e7, 1.0e12          # faixa de validade declarada


def gqrs(E, cc=True):
    return (GQRS_CC if cc else GQRS_NC)*E**GQRS_ALPHA


def carrega(p, col=2):
    d = np.loadtxt(p)
    return d[:, 0], d[:, col]


def carrega_csv(p):
    linhas = [l for l in open(p) if not l.startswith("#") and l.strip()]
    cab = linhas[0].strip().split(",")
    dados = np.array([[float(v) for v in l.split(",")] for l in linhas[1:]])
    return {c: dados[:, i] for i, c in enumerate(cab)}


D = "dipole/data/"
E_cc_g, s_cc_g = carrega(D + "sigma_nuN_CC_GBW.dat")
E_cc_i, s_cc_i = carrega(D + "sigma_nuN_CC_IIM.dat")
E_nc_g, s_nc_g = carrega(D + "sigma_nuN_NC_GBW.dat")
E_nc_i, s_nc_i = carrega(D + "sigma_nuN_NC_IIM.dat")
E_col, s_col = carrega("data/sigma_nuN_CC_collinear.dat")


def faixa(ax):
    ax.axvspan(GQRS_LO, GQRS_HI, color=C["aux"], alpha=0.10, lw=0, zorder=0)


# =====================================================================
# 1. sigma(E), CC e NC, contra a lei de potencia da literatura
# =====================================================================
fig, ax = plt.subplots(2, 2, figsize=(8.6, 6.4), sharex=True,
                       gridspec_kw={"height_ratios": [2.4, 1], "hspace": 0.06,
                                    "wspace": 0.22})

for k, (titulo, dados, lit) in enumerate([
        ("Corrente carregada",
         [("dipolo GBW", E_cc_g, s_cc_g, C["gbw"]),
          ("dipolo bCGC", E_cc_i, s_cc_i, C["iim"]),
          ("colinear NLO (yadism + NNPDF3.1)", E_col, s_col, C["col"])],
         True),
        ("Corrente neutra",
         [("dipolo GBW", E_nc_g, s_nc_g, C["gbw"]),
          ("dipolo bCGC", E_nc_i, s_nc_i, C["iim"])],
         False)]):

    a, b = ax[0, k], ax[1, k]
    faixa(a); faixa(b)

    Eg = np.logspace(np.log10(GQRS_LO), np.log10(GQRS_HI), 200)
    a.plot(Eg, gqrs(Eg, lit), color=C["lit"], ls="--", lw=1.6, zorder=5,
           label=r"Gandhi et al., $%.2f\times10^{-36}\,E^{0{,}363}$"
                 % ((GQRS_CC if lit else GQRS_NC)*1e36))

    for nome, E, s, c in dados:
        a.plot(E, s, color=c, lw=1.5, label=nome)
        m = (E >= GQRS_LO) & (E <= GQRS_HI)
        b.plot(E[m], s[m]/gqrs(E[m], lit), color=c, lw=1.5)

    b.axhline(1.0, color=C["lit"], ls="--", lw=1.0)
    a.set_xscale("log"); a.set_yscale("log")
    a.set_title(titulo, fontsize=10, pad=6)
    a.legend(loc="upper left")
    b.set_xscale("log"); b.set_ylim(0.30, 1.25)
    b.set_xlabel(r"$E_\nu$  [GeV]")
    a.xaxis.set_minor_formatter(NullFormatter())

ax[0, 0].set_ylabel(r"$\sigma_{\nu N}$  [cm$^2$]")
ax[1, 0].set_ylabel("modelo / Gandhi et al.")
ax[0, 0].set_ylim(1e-36, 5e-31)
ax[0, 1].set_ylim(1e-36, 5e-31)
ax[1, 0].annotate("o colinear moderno tambem fica\nabaixo: nem toda a diferenca\ne saturacao",
                  xy=(6e11, 0.63), xytext=(3.0e3, 1.06), fontsize=7.6, color="#555",
                  arrowprops=dict(arrowstyle="->", lw=0.7, color="#888",
                                  connectionstyle="arc3,rad=0.2"))
fig.text(0.5, -0.015,
         "faixa sombreada: onde a lei de potencia e declarada valida "
         r"($10^{16}$ a $10^{21}$ eV). Fora dela a curva nao e desenhada. "
         "O ajuste e de 1996, com PDFs da epoca.",
         ha="center", fontsize=8, color="#555")
fig.savefig("plots/sigma_vs_literatura.pdf")
fig.savefig("plots/sigma_vs_literatura.png")
print("plots/sigma_vs_literatura.{pdf,png}")

# =====================================================================
# 2. razao NC/CC
# =====================================================================
fig, ax = plt.subplots(figsize=(6.4, 4.0))
faixa(ax)
s2w = 0.23122
gVu, gAu = 0.5 - (4/3)*s2w, 0.5
gVd, gAd = -0.5 + (2/3)*s2w, -0.5
S_NC = 2*(gVu**2 + gAu**2) + 2*(gVd**2 + gAd**2)
contagem = (S_NC/4.0)*(91.1876/80.379)**2

ax.plot(E_cc_g, s_nc_g/s_cc_g, color=C["gbw"], lw=1.6, label="dipolo GBW")
ax.plot(E_cc_i, s_nc_i/s_cc_i, color=C["iim"], lw=1.6, label="dipolo bCGC")
ax.hlines(GQRS_NC/GQRS_CC, GQRS_LO, GQRS_HI, color=C["lit"], ls="--", lw=1.6,
          label=r"Gandhi et al., $2{,}31/5{,}53 = %.4f$" % (GQRS_NC/GQRS_CC))
ax.axhline(contagem, color=C["col"], ls=":", lw=1.5,
           label=r"contagem de cargas, $\frac{\Sigma_{NC}}{\Sigma_{CC}}(M_Z/M_W)^2 = %.4f$"
                 % contagem)
ax.set_xscale("log")
ax.set_xlabel(r"$E_\nu$  [GeV]")
ax.set_ylabel(r"$\sigma_{NC}\,/\,\sigma_{CC}$")
ax.set_ylim(0.28, 0.46)
ax.legend(loc="lower right")
ax.annotate("a subida e o $xF_3$ do CC saindo\nde cena: ele e de valencia, vale\n+70% de $\\sigma_{CC}$ em $10^3$ GeV e\n+0,16% em $10^{14}$",
            xy=(2e4, 0.335), xytext=(1.3e6, 0.302), fontsize=8, color="#555",
            arrowprops=dict(arrowstyle="->", lw=0.7, color="#888"))
fig.savefig("plots/razao_nc_cc.pdf"); fig.savefig("plots/razao_nc_cc.png")
print("plots/razao_nc_cc.{pdf,png}")

# =====================================================================
# 3. alpha_eff e o teto 3^(alpha/2)
# =====================================================================
g = carrega_csv("data/alpha_scan_GBW.csv")
i = carrega_csv("data/alpha_scan_IIM.csv")

fig, ax = plt.subplots(2, 1, figsize=(6.6, 6.0), sharex=True,
                       gridspec_kw={"hspace": 0.07})
for a in ax: faixa(a)

ax[0].plot(g["E_GeV"], g["alpha_dipole"], color=C["gbw"], lw=1.6, label="dipolo GBW")
ax[0].plot(i["E_GeV"], i["alpha_dipole"], color=C["iim"], lw=1.6, label="dipolo bCGC")
ax[0].plot(g["E_GeV"], g["alpha_colinear"], color=C["col"], lw=1.6,
           label="colinear NLO")
ax[0].hlines(GQRS_ALPHA, GQRS_LO, GQRS_HI, color=C["lit"], ls="--", lw=1.6,
             label=r"Gandhi et al., $\alpha \simeq 0{,}363$ (constante)")
ax[0].set_ylabel(r"$\alpha_{\rm eff} = d\ln\sigma\,/\,d\ln E$")
ax[0].set_ylim(0.15, 1.15)
ax[0].legend(loc="upper right")

for nome, d, c in [("GBW", g, C["gbw"]), ("bCGC", i, C["iim"])]:
    ax[1].plot(d["E_GeV"], d["teto_dipole"], color=c, lw=1.6)
ax[1].plot(g["E_GeV"], g["teto_colinear"], color=C["col"], lw=1.6)
ax[1].hlines(3.0**(GQRS_ALPHA/2), GQRS_LO, GQRS_HI, color=C["lit"], ls="--", lw=1.6)
ax[1].set_ylabel(r"teto GR $\times$ saturacao:  $3^{\alpha_{\rm eff}/2}$")
ax[1].set_xlabel(r"$E_\nu$  [GeV]")
ax[1].set_xscale("log")
ax[1].set_ylim(1.08, 1.90)
ax[1].annotate("a QUEDA deste teto e a assinatura;\num reajuste de PDF desloca, nao inclina",
               xy=(2e11, 1.145), xytext=(1.5e8, 1.55), fontsize=8.2, color="#444",
               arrowprops=dict(arrowstyle="->", lw=0.8, color="#888"))
ax[1].annotate(r"$3^{0{,}363/2} = 1{,}219$", xy=(3e9, 1.219), xytext=(3e9, 1.30),
               fontsize=8, color="#333", ha="center",
               arrowprops=dict(arrowstyle="-", lw=0.6, color="#888"))
fig.savefig("plots/alpha_eff.pdf"); fig.savefig("plots/alpha_eff.png")
print("plots/alpha_eff.{pdf,png}")

# =====================================================================
# 4. HERA -- dado de verdade, 239 pontos na janela de saturacao
# =====================================================================
h = np.loadtxt(D + "hera_validation.dat")
Q2, x, y, sr, err, mod = h[:, 0], h[:, 1], h[:, 2], h[:, 3], h[:, 4], h[:, 5]

fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.9),
                       gridspec_kw={"wspace": 0.28})

ax[0].errorbar(sr, mod, xerr=err, fmt="o", ms=3.0, lw=0.6, color=C["gbw"],
               alpha=0.75, mec="none")
lim = [0.15, 1.6]
ax[0].plot(lim, lim, color=C["lit"], ls="--", lw=1.2)
ax[0].set_xlim(lim); ax[0].set_ylim(lim)
ax[0].set_xlabel(r"$\sigma_{\rm red}$  medida (H1 + ZEUS)")
ax[0].set_ylabel(r"$\sigma_{\rm red}$  do dipolo, limite EM")
ax[0].set_aspect("equal")
razao = mod/sr
ax[0].set_title(r"$\langle$modelo/dado$\rangle$ = %.3f   (%d pontos)"
                % (razao.mean(), len(sr)), fontsize=9, pad=6)

sc = ax[1].scatter(x, razao, c=np.log10(Q2), s=11, cmap="viridis", lw=0)
ax[1].axhline(1.0, color=C["lit"], ls="--", lw=1.2)
ax[1].set_xscale("log")
ax[1].set_xlabel(r"$x$")
ax[1].set_ylabel("modelo / dado")
ax[1].set_ylim(0.5, 1.7)
cb = fig.colorbar(sc, ax=ax[1], pad=0.02)
cb.set_label(r"$\log_{10}(Q^2/{\rm GeV}^2)$", fontsize=8.5)
fig.text(0.5, -0.06,
         "H1 e ZEUS combinados, Eur. Phys. J. C 75 (2015) 580 [arXiv:1506.06042]. "
         r"Janela $x < 10^{-2}$, $0{,}25 \leq Q^2 \leq 50$ GeV$^2$: "
         "onde a saturacao vive.", ha="center", fontsize=8, color="#555")
fig.savefig("plots/hera.pdf"); fig.savefig("plots/hera.png")
print("plots/hera.{pdf,png}")

print("\nmedia modelo/dado no HERA: %.4f  (rms dos pulls %.2f)"
      % (razao.mean(), np.sqrt((h[:, 8]**2).mean())))
