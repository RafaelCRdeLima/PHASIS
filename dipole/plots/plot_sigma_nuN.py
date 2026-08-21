import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/sigma_nuN_CC_GBW.dat")

Enu = data[:, 0]
sigma_gev2 = data[:, 1]
sigma_cm2 = data[:, 2]

plt.figure(figsize=(8, 5))

plt.loglog(
    Enu,
    sigma_cm2,
    linewidth=2,
    label=r"$\sigma_{\nu N}^{CC}$ GBW"
)

plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
plt.ylabel(r"$\sigma_{\nu N}^{CC}\ [\mathrm{cm}^2]$")

plt.ylim(1e-38, 1e-28)
plt.xlim(1e3, 1e13)

plt.title(r"Charged-Current Neutrino--Nucleon Cross Section")

# plt.grid(True, which="both", alpha=0.3)

plt.legend()
plt.tight_layout()

plt.savefig("plots/sigma_nuN_CC_GBW.png", dpi=200)

plt.show()