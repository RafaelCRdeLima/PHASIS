import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/integrand.dat")

r = data[:,0]
z = data[:,1]
IT = data[:,2]
IL = data[:,3]

plt.figure(figsize=(8,5))
sc = plt.scatter(r, z, c=np.log10(IT + 1e-300), s=8)
plt.xscale("log")
plt.xlabel(r"$r\ [\mathrm{GeV}^{-1}]$")
plt.ylabel(r"$z$")
plt.colorbar(sc, label=r"$\log_{10}\,\mathcal{I}_T$")
plt.title(r"Integrando transversal: $2\pi r |\Psi_T|^2\sigma_{dip}$")
plt.tight_layout()
plt.savefig("plots/integrandT.png", dpi=200)
plt.show()

plt.figure(figsize=(8,5))
sc = plt.scatter(r, z, c=np.log10(IL + 1e-300), s=8)
plt.xscale("log")
plt.xlabel(r"$r\ [\mathrm{GeV}^{-1}]$")
plt.ylabel(r"$z$")
plt.colorbar(sc, label=r"$\log_{10}\,\mathcal{I}_L$")
plt.title(r"Integrando longitudinal: $2\pi r |\Psi_L|^2\sigma_{dip}$")
plt.tight_layout()
plt.savefig("plots/integrandL.png", dpi=200)
plt.show()