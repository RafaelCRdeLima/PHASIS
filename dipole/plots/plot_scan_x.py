import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/scan_x_F2_CC.dat")

x  = data[:, 0]
Q2 = data[:, 1]
FT = data[:, 2]
FL = data[:, 3]
F2 = data[:, 4]

plt.figure(figsize=(8, 5))

plt.plot(x, FT, label=r"$F_T$")
plt.plot(x, FL, label=r"$F_L$")
plt.plot(x, F2, label=r"$F_2 = F_T + F_L$", linewidth=2)

plt.xscale("log")

plt.xlabel(r"$x$")
plt.ylabel(r"$F_i(x,Q^2)$")
plt.title(fr"Scan em $x$ para $Q^2={Q2[0]:.2e}\,\mathrm{{GeV}}^2$")

plt.legend()
#plt.grid(True, which="both", alpha=0.3)
plt.xlim(1e-8,7e-3)
plt.tight_layout()

plt.savefig("plots/scan_x_F2_CC.png", dpi=200)
plt.show()