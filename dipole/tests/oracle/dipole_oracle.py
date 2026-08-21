"""Oraculo Python do calculo de sigma_nuN CC.

Reimplementacao fiel de src/{wavefunctions,dipole_models,integrals,sigma_nuN}.cpp,
incluindo as grades de quadratura originais.

PROPOSITO: este modulo reproduz o comportamento *atual* do C++, nao o correto.
Serve para provar que uma mudanca no C++ fez exatamente o que se pretendia.
Concordancia medida com build/sigma_nuN: 3e-4 (sem o termo F3, que vale ~3e-4
em 1e6 GeV).

Nao inclui xF3 (exigiria LHAPDF); rode o C++ com --use-F3 0 para comparacao exata.
"""

import numpy as np
from scipy.special import kv

pi = np.pi

# ---- constantes: parameters.hpp -------------------------------------------
GF = 1.1663787e-5      # GeV^-2
MW = 80.379            # GeV
MZ = 91.1876           # GeV
MN = 0.938272          # GeV
GEV2_TO_CM2 = 0.389379e-27

# CC: gV=-1, gA=1, alphaEW = (MW^2 GF/sqrt2)/(4pi)   [parameters.hpp:81-85]
ALPHA_CC = (MW * MW * GF / np.sqrt(2.0)) / (4.0 * pi)

# ---- GBW: dipole_models.hpp:6-9 -------------------------------------------
GBW_SIGMA0_MB = 29.12
GBW_LAMBDA = 0.277
GBW_X0 = 0.41e-4
GBW_Q0SQ = 1.0
GBW_SIGMA0 = GBW_SIGMA0_MB / 0.389379   # GeV^-2

# ---- massas de quark: parameters.hpp:19-26 --------------------------------
M_U = M_D = M_S = 0.03
M_C = 1.3113

# canais somados por computeFCCtotal: (u,d) e (c,s)
CC_CHANNELS = ((M_U, M_D), (M_C, M_S))


def sigma_dip_gbw(r, x):
    """dipole_models.cpp:R0sq_GBW + sigmaDipoleGBW. Sem corte de validade em x."""
    R0sq = (1.0 / GBW_Q0SQ) * (x / GBW_X0) ** GBW_LAMBDA
    return GBW_SIGMA0 * (1.0 - np.exp(-r * r / (4.0 * R0sq)))


def eps2(z, Q2, m, mu):
    """wavefunctions.cpp:epsilon2"""
    return z * (1 - z) * Q2 + z * m * m + (1 - z) * mu * mu


def psiT2(r, z, Q2, m, mu, gV2=1.0, gA2=1.0, alpha=ALPHA_CC):
    """wavefunctions.cpp:psiT2"""
    e2 = eps2(z, Q2, m, mu)
    a = np.sqrt(e2) * r
    K0, K1 = kv(0, a), kv(1, a)
    zb = 1 - z
    termK1 = (gV2 + gA2) * (z * z + zb * zb) * e2 * K1 * K1
    termK0 = (gV2 * (z * m + zb * mu) ** 2 + gA2 * (z * m - zb * mu) ** 2) * K0 * K0
    return 6.0 * alpha / (2.0 * pi) ** 2 * (termK1 + termK0)


def psiL2(r, z, Q2, m, mu, gV2=1.0, gA2=1.0, alpha=ALPHA_CC):
    """wavefunctions.cpp:psiL2"""
    e2 = eps2(z, Q2, m, mu)
    a = np.sqrt(e2) * r
    K0, K1 = kv(0, a), kv(1, a)
    zb = 1 - z
    A = gV2 * (m - mu) ** 2 + gA2 * (m + mu) ** 2
    Bv = 2 * Q2 * z * zb + (m - mu) * (z * m - zb * mu)
    Ba = 2 * Q2 * z * zb + (m + mu) * (z * m + zb * mu)
    return 6.0 * alpha / ((2.0 * pi) ** 2 * Q2) * (
        A * e2 * K1 * K1 + (gV2 * Bv * Bv + gA2 * Ba * Ba) * K0 * K0
    )


def simpson_nodes(a, b, N):
    """integrals.cpp:simpson - nos e pesos."""
    if N % 2:
        N += 1
    h = (b - a) / N
    xs = a + h * np.arange(N + 1)
    w = np.ones(N + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return xs, w * h / 3.0


def F_TL(x, Q2, m, mu, Nr, Nz, rmin=1e-6, rmax=1e2, zmin=1e-6, zmax=1 - 1e-6,
         logr=False, logz=False):
    """integrals.cpp:FT_GBW/FL_GBW, ja divididos por alphaEW (sigma_nuN.cpp:75).

    logr/logz=False reproduz o codigo original (grades uniformes).
    """
    if logr:
        u, wu = simpson_nodes(np.log(rmin), np.log(rmax), Nr)
        r = np.exp(u)
        wr = wu * r
    else:
        r, wr = simpson_nodes(rmin, rmax, Nr)
    if logz:
        half = max(2, Nz // 2)
        t, wt = simpson_nodes(np.log(zmin), np.log(0.5), half)
        zl = np.exp(t)
        z = np.concatenate([zl, 1.0 - zl[::-1]])
        wz = np.concatenate([wt * zl, (wt * zl)[::-1]])
    else:
        z, wz = simpson_nodes(zmin, zmax, Nz)

    R, Z = r[:, None], z[None, :]
    sd = sigma_dip_gbw(R, x)
    iT = 2 * pi * R * psiT2(R, Z, Q2, m, mu) * sd
    iL = 2 * pi * R * psiL2(R, Z, Q2, m, mu) * sd
    W = wr[:, None] * wz[None, :]
    c = Q2 / (4 * pi * pi) / ALPHA_CC
    return c * (W * iT).sum(), c * (W * iL).sum()


def large_x_factor(x):
    """sigma_nuN.cpp:largeXFactor"""
    return (1.0 - x) ** 7.0 if 0.0 < x < 1.0 else 0.0


def F2_total(x, Q2, Nr, Nz, **kw):
    """sigma_nuN.cpp:computeFCCtotal - soma canais (u,d) + (c,s), aplica (1-x)^7."""
    FT = FL = 0.0
    for m, mu in CC_CHANNELS:
        t, l = F_TL(x, Q2, m, mu, Nr, Nz, **kw)
        FT += t
        FL += l
    lx = large_x_factor(x)
    return lx * FT, lx * FL


def sigma_nuN_CC(E, NlogQ=16, Nlogx=16, Nr=50, Nz=30, Q2min=1.0, **kw):
    """sigma_nuN.cpp:sigmaNuN_CC. Defaults = os de main(), que o Makefile usa.

    Retorna GeV^-2. Sem o termo xF3.
    """
    s = 2 * MN * E
    Q2max = 0.999 * s
    if Q2max <= Q2min:
        return 0.0
    lq, wq = simpson_nodes(np.log(Q2min), np.log(Q2max), NlogQ)
    total = 0.0
    for lQ, wQ in zip(lq, wq):
        Q2 = np.exp(lQ)
        xmin = Q2 / s
        xmax = 0.999999
        if xmin >= xmax:
            continue
        lxs, wx = simpson_nodes(np.log(xmin), np.log(xmax), Nlogx)
        inner = 0.0
        for lx, w in zip(lxs, wx):
            x = np.exp(lx)
            y = Q2 / (x * s)
            if not (0.0 < y < 1.0):
                continue
            FT, FL = F2_total(x, Q2, Nr, Nz, **kw)
            F2 = FT + FL
            prop = (MW * MW / (Q2 + MW * MW)) ** 2
            pre = GF * GF * MN * E / pi * prop
            bracket = 0.5 * (1 + (1 - y) ** 2) * F2 - 0.5 * y * y * FL
            inner += w * (pre * bracket) / s
        total += wQ * Q2 * inner
    return total


def sigma_cm2(E, **kw):
    return sigma_nuN_CC(E, **kw) * GEV2_TO_CM2
