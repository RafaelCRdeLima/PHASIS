import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------
# leitura do arquivo
# ---------------------------------

data = np.loadtxt("data/wavefunctions.dat")

r     = data[:,0]
z     = data[:,1]
psiT2 = data[:,2]
psiL2 = data[:,3]

# ---------------------------------
# mapa de calor transversal
# ---------------------------------

plt.figure(figsize=(8,5))

sc = plt.scatter(
    r,
    z,
    c=np.log10(psiT2 + 1e-300),
    s=8
)

plt.xscale("log")

plt.xlabel(r"$r\ [\mathrm{GeV}^{-1}]$")
plt.ylabel(r"$z$")

cbar = plt.colorbar(sc)
cbar.set_label(r"$\log_{10} |\Psi_T|^2$")

plt.title(r"Transverse wavefunction")

plt.tight_layout()

plt.savefig("plots/psiT2_heatmap.png", dpi=200)

plt.show()

# ---------------------------------
# mapa de calor longitudinal
# ---------------------------------

plt.figure(figsize=(8,5))

sc = plt.scatter(
    r,
    z,
    c=np.log10(psiL2 + 1e-300),
    s=8
)

plt.xscale("log")

plt.xlabel(r"$r\ [\mathrm{GeV}^{-1}]$")
plt.ylabel(r"$z$")

cbar = plt.colorbar(sc)
cbar.set_label(r"$\log_{10} |\Psi_L|^2$")

plt.title(r"Longitudinal wavefunction")

plt.tight_layout()

plt.savefig("plots/psiL2_heatmap.png", dpi=200)

plt.show()

# =========================================================
# curvas 2D para diferentes z
# =========================================================

# seleciona alguns valores bem separados de z
z_values = np.unique(z)

z_selected = [
    z_values[np.argmin(np.abs(z_values - 0.1))],
    z_values[np.argmin(np.abs(z_values - 0.3))],
    z_values[np.argmin(np.abs(z_values - 0.5))],
    z_values[np.argmin(np.abs(z_values - 0.8))]
]

# ---------------------------------
# curvas transversais
# ---------------------------------

plt.figure(figsize=(8,5))

for z0 in z_selected:

    mask = np.isclose(z, z0)

    r_slice    = r[mask]
    psiT_slice = psiT2[mask]

    order = np.argsort(r_slice)

    plt.plot(
        r_slice[order],
        psiT_slice[order],
        linewidth=2,
        label=rf"$z = {z0:.2f}$"
    )

plt.xscale("log")
plt.yscale("log")

plt.xlabel(r"$r\ [\mathrm{GeV}^{-1}]$")
plt.ylabel(r"$|\Psi_T|^2$")

plt.title(r"Transverse wavefunction for selected $z$")

plt.legend()

plt.tight_layout()

plt.savefig("plots/psiT2_curves.png", dpi=200)

plt.show()

# ---------------------------------
# curvas longitudinais
# ---------------------------------

plt.figure(figsize=(8,5))

for z0 in z_selected:

    mask = np.isclose(z, z0)

    r_slice    = r[mask]
    psiL_slice = psiL2[mask]

    order = np.argsort(r_slice)

    plt.plot(
        r_slice[order],
        psiL_slice[order],
        linewidth=2,
        label=rf"$z = {z0:.2f}$"
    )

plt.xscale("log")
plt.yscale("log")

plt.xlabel(r"$r\ [\mathrm{GeV}^{-1}]$")
plt.ylabel(r"$|\Psi_L|^2$")

plt.title(r"Longitudinal wavefunction for selected $z$")

plt.legend()

plt.tight_layout()

plt.savefig("plots/psiL2_curves.png", dpi=200)

plt.show()