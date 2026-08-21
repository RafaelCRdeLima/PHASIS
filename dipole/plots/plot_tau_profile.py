import sys
import numpy as np
import matplotlib.pyplot as plt


# =====================================================
# Command-line arguments
# =====================================================

input_file = "data/tau_profile.dat"
output_tag = "default"
profile_name = "default"

if len(sys.argv) > 1:
    input_file = sys.argv[1]

if len(sys.argv) > 2:
    output_tag = sys.argv[2]

if len(sys.argv) > 3:
    profile_name = sys.argv[3]
else:
    profile_name = output_tag

profile_title = profile_name.replace("_", " ")
output_tag = output_tag.lower()


# =====================================================
# Load data
# =====================================================

data = np.loadtxt(input_file)

r_cm = data[:, 0]
rho = data[:, 1]
tau = data[:, 2]
survival = data[:, 3]

r_km = r_cm / 1.0e5

tau_nu = 2.0 / 3.0


# =====================================================
# Find neutrinosphere radius: tau(r_nu) = 2/3
# Assumes tau decreases outward.
# =====================================================

idx = np.where((tau[:-1] > tau_nu) & (tau[1:] <= tau_nu))[0]

r_neutrino_km = None

if len(idx) > 0:
    i = idx[0]

    r1, r2 = r_km[i], r_km[i + 1]
    t1, t2 = tau[i], tau[i + 1]

    if t2 != t1:
        r_neutrino_km = r1 + (tau_nu - t1) * (r2 - r1) / (t2 - t1)
    else:
        r_neutrino_km = r1


# =====================================================
# 1) Optical depth profile
# =====================================================

plt.figure(figsize=(8, 5))

plt.semilogy(
    r_km,
    tau,
    label=r"$\tau_\nu(r)$"
)

plt.axhline(
    tau_nu,
    linestyle="--",
    label=r"$\tau_\nu=2/3$"
)

if r_neutrino_km is not None:
    plt.axvline(
        r_neutrino_km,
        linestyle="--",
        label=fr"$r_\nu={r_neutrino_km:.2f}\,\mathrm{{km}}$"
    )

plt.xlabel(r"$r\ [\mathrm{km}]$")
plt.ylabel(r"$\tau_\nu(r)$")
plt.title(fr"Optical Depth Profile — {profile_title}")

plt.legend(fontsize=9)
plt.tight_layout()

plt.savefig(
    f"plots/tau_profile_neutrinosphere_{output_tag}.png",
    dpi=200
)

plt.show()


# =====================================================
# 2) Survival probability profile
# =====================================================

plt.figure(figsize=(8, 5))

plt.plot(
    r_km,
    survival,
    label=r"$e^{-\tau_\nu(r)}$"
)

if r_neutrino_km is not None:
    plt.axvline(
        r_neutrino_km,
        linestyle="--",
        label=fr"$r_\nu={r_neutrino_km:.2f}\,\mathrm{{km}}$"
    )

plt.xlabel(r"$r\ [\mathrm{km}]$")
plt.ylabel(r"Survival probability")
plt.title(fr"Neutrino Escape Region — {profile_title}")

plt.legend(fontsize=9)
plt.tight_layout()

plt.savefig(
    f"plots/survival_profile_neutrinosphere_{output_tag}.png",
    dpi=200
)

plt.show()


# =====================================================
# Console summary
# =====================================================

if r_neutrino_km is not None:
    print(f"Neutrinosphere radius: r_nu = {r_neutrino_km:.6e} km")
else:
    print("No tau = 2/3 crossing found in this profile.")