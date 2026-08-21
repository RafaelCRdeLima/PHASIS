import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/scan_x_F2_CC_models.dat")

x = data[:, 0]
Q2 = data[:, 1]

FT_GBW = data[:, 2]
FL_GBW = data[:, 3]
F2_GBW = data[:, 4]

xF3 = data[:, 5]

FT_IIM = data[:, 6]
FL_IIM = data[:, 7]
F2_IIM = data[:, 8]

plt.figure(figsize=(8, 5))

plt.plot(x, F2_GBW, label=r"$F_2$ GBW", linewidth=2)
plt.plot(x, F2_IIM, label=r"$F_2$ IIM", linewidth=2)
plt.plot(x, xF3, label=r"$xF_3$", linewidth=2, linestyle="--")

plt.xscale("log")
plt.yscale("log")

plt.xlabel(r"$x$")
plt.ylabel(r"Structure Functions")
plt.title(
    fr"GBW vs IIM with $xF_3$, "
    fr"$Q^2={Q2[0]:.2e}\,\mathrm{{GeV}}^2$"
)

#plt.grid(True, which="both", alpha=0.3)
plt.legend()
plt.tight_layout()
plt.xlim(1e-8, 7e-3)

plt.savefig("plots/scan_x_F2_CC_models.png", dpi=200)
plt.show()