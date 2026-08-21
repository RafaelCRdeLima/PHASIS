import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterSciNotation
from matplotlib.ticker import FuncFormatter


DATA_DIR = "data"
PLOT_DIR = "plots"
os.makedirs(PLOT_DIR, exist_ok=True)


def sci_label(value):
    mantissa, exponent = f"{value:.0e}".split("e")
    return rf"${mantissa}\times10^{{{int(exponent)}}}$"


# ============================================================
# 1. Neutrino--Nucleon Cross Section: GBW vs IIM
# ============================================================

sigma_gbw = np.loadtxt(f"{DATA_DIR}/sigma_nuN_CC_GBW.dat")
sigma_iim = np.loadtxt(f"{DATA_DIR}/sigma_nuN_CC_IIM.dat")

E_gbw = sigma_gbw[:, 0]
sigma_gbw_cm2 = sigma_gbw[:, 1]

E_iim = sigma_iim[:, 0]
sigma_iim_cm2 = sigma_iim[:, 1]

plt.figure(figsize=(8, 5))

plt.loglog(E_gbw, sigma_gbw_cm2, label="GBW", linewidth=2)
plt.loglog(E_iim, sigma_iim_cm2, label="IIM", linewidth=2)

plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
plt.ylabel(r"$\sigma_{\nu N}^{\mathrm{CC}}\ [\mathrm{cm}^2]$")
plt.title("Charged-Current Neutrino--Nucleon Cross Section")
plt.legend()
plt.tight_layout()

plt.savefig(f"{PLOT_DIR}/sigma_nuN_CC_GBW_vs_IIM.png", dpi=300)
plt.show()


# ============================================================
# 2. Optical Depth Profiles: GBW vs IIM
# ============================================================

profiles_gbw = np.loadtxt(f"{DATA_DIR}/tau_profiles_scan_GBW.dat")
profiles_iim = np.loadtxt(f"{DATA_DIR}/tau_profiles_scan_IIM.dat")

selected = [1e5, 1e7, 1e9, 1e12]

plt.figure(figsize=(8, 5))

for E in selected:
    mask_gbw = np.isclose(profiles_gbw[:, 0], E, rtol=1e-10, atol=0.0)
    mask_iim = np.isclose(profiles_iim[:, 0], E, rtol=1e-10, atol=0.0)

    plt.semilogy(
        profiles_gbw[mask_gbw, 2],
        profiles_gbw[mask_gbw, 3],
        linewidth=2,
        label=rf"GBW, $E_\nu={sci_label(E)}\,\mathrm{{GeV}}$"
    )

    plt.semilogy(
        profiles_iim[mask_iim, 2],
        profiles_iim[mask_iim, 3],
        linewidth=2,
        linestyle="--",
        label=rf"IIM, $E_\nu={sci_label(E)}\,\mathrm{{GeV}}$"
    )

plt.axhline(2.0/3.0, linestyle=":", color="black", label=r"$\tau=2/3$")

plt.xlabel(r"$r\ [\mathrm{km}]$")
plt.ylabel(r"$\tau_\nu(r,E_\nu)$")
plt.title("Optical Depth Profiles")
plt.legend(fontsize=8)
plt.tight_layout()

energy_tag = "_".join([f"E{E:.0e}".replace("+", "") for E in selected])
plt.savefig(f"{PLOT_DIR}/tau_profiles_GBW_vs_IIM_{energy_tag}.png", dpi=300)
plt.show()


# ============================================================
# 3. UHE Neutrinosphere Radius: GBW vs IIM
# ============================================================

sphere_gbw = np.loadtxt(f"{DATA_DIR}/neutrinosphere_scan_GBW.dat")
sphere_iim = np.loadtxt(f"{DATA_DIR}/neutrinosphere_scan_IIM.dat")

E_gbw = sphere_gbw[:, 0]
R_gbw_km = sphere_gbw[:, 2]

E_iim = sphere_iim[:, 0]
R_iim_km = sphere_iim[:, 2]

plt.figure(figsize=(8, 5))

scale = 1e9

plt.plot(E_gbw, R_gbw_km / scale, label="GBW", linewidth=2)
plt.plot(E_iim, R_iim_km / scale, label="IIM", linewidth=2)

plt.xscale("log")

plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
plt.ylabel(r"$R_{\tau=1}\ [10^9\,\mathrm{km}]$")
plt.title("UHE Neutrinosphere Radius")

plt.gca().yaxis.set_major_formatter(
    FuncFormatter(lambda y, _: f"{y:g}")
)

plt.legend()
plt.tight_layout()

plt.savefig(
    f"{PLOT_DIR}/neutrinosphere_radius_GBW_vs_IIM.png",
    dpi=300
)

plt.show()
