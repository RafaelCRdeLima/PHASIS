import sys
import numpy as np
import matplotlib.pyplot as plt


# =====================================================
# Command-line arguments
# =====================================================

model = "GBW"
profile_name = "default"

if len(sys.argv) > 1:
    model = sys.argv[1]

if len(sys.argv) > 2:
    profile_name = sys.argv[2]

model_tag = model.lower()
profile_tag = profile_name.lower()
profile_title = profile_name.replace("_", " ")


# =====================================================
# Input files
# =====================================================

summary_file = f"data/neutrinosphere_scan_{model}.dat"
profiles_file = f"data/tau_profiles_scan_{model}.dat"


# =====================================================
# Load scan summary
# =====================================================

summary = np.loadtxt(summary_file)

Enu = summary[:, 0]
sigma = summary[:, 1]
rnu_km = summary[:, 3]
tau_center = summary[:, 4]


# =====================================================
# 1) sigma_nuN(E_nu)
# =====================================================

plt.figure(figsize=(8, 5))
plt.loglog(Enu, sigma, marker="o", label=model)

plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
plt.ylabel(r"$\sigma_{\nu N}^{CC}\ [\mathrm{cm}^{2}]$")
plt.title(
    fr"Neutrino-Nucleon Cross Section — {model}"
)

plt.legend()
plt.tight_layout()

plt.savefig(
    f"plots/scan_sigma_nuN_{model_tag}_{profile_tag}.png",
    dpi=200
)

plt.show()


# =====================================================
# 2) r_nu(E_nu)
# =====================================================

plt.figure(figsize=(8, 5))
plt.semilogx(Enu, rnu_km, marker="o", label=model)

plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
plt.ylabel(r"$r_\nu\ [\mathrm{km}]$")
plt.title(
    fr"UHE Neutrinosphere Radius — {model} ({profile_title})"
)

plt.legend()
plt.tight_layout()

plt.savefig(
    f"plots/scan_rnu_{model_tag}_{profile_tag}.png",
    dpi=200
)

plt.show()


# =====================================================
# 3) tau(r,E_nu) for selected energies
# =====================================================

profiles = np.loadtxt(profiles_file)

E_prof = profiles[:, 0]
r_km = profiles[:, 2]
tau = profiles[:, 3]

selected = [1e5, 1e7, 1e9, 1e12]

energy_tag = "_".join([
    f"E{E:.0e}".replace("+", "")
    for E in selected
])

plt.figure(figsize=(8, 5))

for E in selected:
    mask = np.isclose(E_prof, E, rtol=1e-6, atol=0.0)

    if not np.any(mask):
        print(f"Warning: energy E = {E:.3e} GeV not found in {profiles_file}")
        continue

    mantissa, exponent = f"{E:.0e}".split("e")

    label = (
        fr"$E_\nu={mantissa}\times10^{{{int(exponent)}}}"
        fr"\,\mathrm{{GeV}}$"
    )

    plt.semilogy(
        r_km[mask],
        tau[mask],
        label=label
    )

plt.axhline(
    2.0/3.0,
    linestyle="--",
    label=r"$\tau=2/3$"
)

plt.xlabel(r"$r\ [\mathrm{km}]$")
plt.ylabel(r"$\tau_\nu(r,E_\nu)$")
plt.title(
    fr"Radial Optical Depth Profiles — {model} ({profile_title})"
)

plt.legend(fontsize=9)
plt.tight_layout()

plt.savefig(
    f"plots/scan_tau_profiles_{model_tag}_{profile_tag}_{energy_tag}.png",
    dpi=200
)

plt.show()