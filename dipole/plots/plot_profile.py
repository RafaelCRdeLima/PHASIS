import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/flash_like_profile.dat")

r_cm = data[:, 0]
dens = data[:, 1]
temp = data[:, 2]
ye   = data[:, 3]
pres = data[:, 4]
entr = data[:, 5]

r_km = r_cm / 1.0e5

plt.figure(figsize=(8, 5))
plt.semilogy(r_km, dens)
plt.xlabel(r"$r\ [\mathrm{km}]$")
plt.ylabel(r"$\rho\ [\mathrm{g\,cm^{-3}}]$")
plt.title("Perfil de densidade")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("plots/profile_density.png", dpi=200)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(r_km, temp)
plt.xlabel(r"$r\ [\mathrm{km}]$")
plt.ylabel(r"$T\ [\mathrm{MeV}]$")
plt.title("Perfil de temperatura")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("plots/profile_temperature.png", dpi=200)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(r_km, ye)
plt.xlabel(r"$r\ [\mathrm{km}]$")
plt.ylabel(r"$Y_e$")
plt.title("Perfil de fração eletrônica")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("plots/profile_ye.png", dpi=200)
plt.show()

plt.figure(figsize=(8, 5))
plt.semilogy(r_km, pres)
plt.xlabel(r"$r\ [\mathrm{km}]$")
plt.ylabel(r"$P_{\rm dummy}$")
plt.title("Perfil de pressão dummy")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("plots/profile_pressure.png", dpi=200)
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(r_km, entr)
plt.xlabel(r"$r\ [\mathrm{km}]$")
plt.ylabel(r"$s_{\rm dummy}$")
plt.title("Perfil de entropia dummy")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("plots/profile_entropy.png", dpi=200)
plt.show()