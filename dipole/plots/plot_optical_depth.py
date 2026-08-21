import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/optical_depth_CC_GBW.dat")

Enu = data[:, 0]
sigma = data[:, 1]
tau = data[:, 2]
survival = data[:, 3]

plt.figure(figsize=(8, 5))
plt.loglog(Enu, tau, marker="o", label=r"$\tau_\nu(E_\nu)$")
plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
plt.ylabel(r"$\tau_\nu$")
plt.title(r"Profundidade óptica neutrino-matéria")
plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("plots/optical_depth_CC_GBW.png", dpi=200)
plt.show()

log10_survival = -tau / np.log(10)

plt.figure(figsize=(8, 5))
plt.semilogx(Enu, log10_survival, marker="o")
plt.xlabel(r"$E_\nu\ [\mathrm{GeV}]$")
plt.ylabel(r"$\log_{10}(e^{-\tau_\nu})$")
plt.title(r"Log da probabilidade de sobrevivência")
plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("plots/log_survival_probability_CC_GBW.png", dpi=200)
plt.show()