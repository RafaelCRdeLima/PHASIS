"""sigma_CC(E) colinear, pela MESMA formula mestra que o dipolo usa.

A comparacao dipolo x colinear so mede as FUNCOES DE ESTRUTURA se tudo
o mais for identico. Entao esta rotina reproduz linha por linha
dipole/src/sigma_nuN_core.cpp:

  sigma = INT dQ^2 INT dx  (1/(x s)) d2sigma/dx dy
  d2sigma/dx dy = (G_F^2 M_N E / pi) (M_W^2/(Q^2+M_W^2))^2
                  [ (1+(1-y)^2)/2 F2 - y^2/2 F_L + s3 y(1-y/2) xF3 ]

com y = Q^2/(x s), s3 = +1 para nu, quadratura de Simpson composta em
ln Q^2 e ln x, e a MESMA regra de densidade de nos por decada
(nodesForRange), inclusive o Nlogx fixo ao longo do laco de Q^2.

DUAS diferencas deliberadas, e so duas:

 1. As funcoes de estrutura sao colineares (yadism NLO + NNPDF3.1),
    nao do modelo de dipolo. E o que se quer medir.

 2. Sem o fator (1-x)^7. Ele e a regra de contagem de constituintes de
    Kutak-Kwiecinski, uma prescricao do modelo de dipolo para grande x;
    os PDFs colineares ja se anulam em x -> 1 por conta propria, e
    aplica-lo aqui seria contar duas vezes.

y = 1 e tratado como EXTREMO CINEMATICO, com valor finito, igual ao C++
depois do conserto -- zera-lo punha um degrau no extremo da quadratura e
custava um erro O(h) em vez de O(h^4).
"""
import argparse
import subprocess

import numpy as np
from scipy.interpolate import RectBivariateSpline

import theory_card as TC


def nodes_for_range(log_lo, log_hi, per_decade, n_min=16):
    """dipole/src/sigma_nuN_core.cpp:nodesForRange"""
    decades = (log_hi - log_lo)/np.log(10.0)
    n = int(np.ceil(per_decade*decades))
    if n < n_min:
        n = n_min
    if n % 2:
        n += 1
    return n


def simpson_nodes(a, b, n):
    """dipole/src/integrals.cpp:simpsonRule"""
    if n < 2:
        n = 2
    if n % 2:
        n += 1
    h = (b - a)/n
    xs = a + h*np.arange(n + 1)
    w = np.ones(n + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return xs, w*h/3.0


class SFTable:
    """F2, FL, xF3 interpolados em (ln x, ln Q^2), com truncamento nas
    bordas -- sem extrapolar, igual a StructureTable do dipole."""

    def __init__(self, npz):
        d = np.load(npz)
        self.lx = np.log(d["x"])
        self.lq = np.log(d["q2"])
        self.x_min, self.x_max = d["x"][0], d["x"][-1]
        self.q2_min, self.q2_max = d["q2"][0], d["q2"][-1]
        self.pdf_set = str(d["pdf_set"])
        self.pto = int(d["pto"])
        self.yadism = str(d["yadism_version"])
        self.pdf_xmin = float(d["pdf_xmin"]) if "pdf_xmin" in d else TC.PDF_XMIN
        # spline bicubica em (ln Q^2, ln x) sobre os valores CRUS: F_L
        # pode ser negativa em NLO perto de x -> 1, e log dela nao existe
        self._f2 = RectBivariateSpline(self.lq, self.lx, d["F2"], kx=3, ky=3, s=0)
        self._fl = RectBivariateSpline(self.lq, self.lx, d["FL"], kx=3, ky=3, s=0)
        self._f3 = RectBivariateSpline(self.lq, self.lx, d["xF3"], kx=3, ky=3, s=0)

    def at(self, x, q2):
        lx = np.clip(np.log(np.asarray(x, dtype=float)), self.lx[0], self.lx[-1])
        lq = np.clip(np.log(np.asarray(q2, dtype=float)), self.lq[0], self.lq[-1])
        lq = np.broadcast_to(lq, lx.shape)
        return (self._f2(lq, lx, grid=False),
                self._fl(lq, lx, grid=False),
                self._f3(lq, lx, grid=False))


def d2sigma_dxdy(E, x, q2, sf, s3=+1.0):
    s = 2.0*TC.MN*E
    y = q2/(x*s)
    out = np.zeros_like(x, dtype=float)

    ok = (x > 0.0) & (x < 1.0) & (y > 0.0)
    # y = 1 e o extremo cinematico, nao um ponto proibido. Acima de
    # 1 + 1e-9 e fora da cinematica de verdade; ate la e so o
    # arredondamento de exp(log(Q^2/s)) no no de borda.
    ok &= (y <= 1.0 + 1.0e-9)
    if not ok.any():
        return out
    yy = np.minimum(y[ok], 1.0)
    xx = x[ok]

    F2, FL, xF3 = sf.at(xx, q2)
    prop = (TC.MW**2/(q2 + TC.MW**2))**2
    pre = TC.GF**2*TC.MN*E/np.pi*prop
    bracket = (0.5*(1.0 + (1.0 - yy)**2)*F2
               - 0.5*yy*yy*FL
               + s3*yy*(1.0 - 0.5*yy)*xF3)
    out[ok] = pre*bracket
    return out


def sigma_cc(E, sf, per_dec_q=16.0, per_dec_x=16.0, q2min=1.0,
             conta_extrapolado=False):
    """dipole/src/sigma_nuN_core.cpp:sigmaNuN_CC, em GeV^-2."""
    s = 2.0*TC.MN*E
    q2max = 0.999*s
    if q2max <= q2min:
        return (0.0, 0.0) if conta_extrapolado else 0.0

    lq, wq = simpson_nodes(np.log(q2min), np.log(q2max),
                           nodes_for_range(np.log(q2min), np.log(q2max), per_dec_q))
    xmax = 0.999999
    # FIXO ao longo do laco de Q^2, pela faixa mais larga -- ver o
    # comentario em sigma_nuN_core.cpp sobre o ruido que o contrario
    # injetava.
    n_x = nodes_for_range(np.log(q2min/s), np.log(xmax), per_dec_x)

    total = 0.0
    extrap = 0.0
    for lQ, wQ in zip(lq, wq):
        q2 = np.exp(lQ)
        xmin = q2/s
        if xmin >= xmax:
            continue
        lxs, wx = simpson_nodes(np.log(xmin), np.log(xmax), n_x)
        xv = np.exp(lxs)
        g = d2sigma_dxdy(E, xv, q2, sf)/s
        contrib = wQ*q2*np.sum(wx*g)
        total += contrib
        if conta_extrapolado:
            m = xv < sf.pdf_xmin
            if m.any():
                extrap += wQ*q2*np.sum(wx[m]*g[m])

    return (total, extrap/total if total else 0.0) if conta_extrapolado else total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sf", default="data/collinear_sf.npz")
    ap.add_argument("--out", default="data/sigma_nuN_CC_collinear.dat")
    ap.add_argument("--logEmin", type=float, default=3.0)
    ap.add_argument("--logEmax", type=float, default=14.0)
    ap.add_argument("--NE", type=int, default=300)
    ap.add_argument("--nodesQ", type=float, default=16.0)
    ap.add_argument("--nodesX", type=float, default=16.0)
    ap.add_argument("--q2min", type=float, default=1.0)
    a = ap.parse_args()

    sf = SFTable(a.sf)
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    except Exception:
        commit = "desconhecido"

    Es = np.logspace(a.logEmin, a.logEmax, a.NE)
    sig = np.empty(a.NE)
    frac = np.empty(a.NE)
    for i, E in enumerate(Es):
        sig[i], frac[i] = sigma_cc(E, sf, a.nodesQ, a.nodesX, a.q2min,
                                   conta_extrapolado=True)
        if i % 25 == 0:
            print(f"  E = {E:10.4g} GeV   sigma = {sig[i]*TC.GEV2_TO_CM2:.5e} cm^2"
                  f"   x<{sf.pdf_xmin:g}: {100*frac[i]:.1f}% de sigma", flush=True)

    with open(a.out, "w") as f:
        f.write("# convention_y = (E_in - E_out)/E_in\n")
        f.write("# target       = isoscalar_nucleon\n")
        f.write("# projectile   = nu\n")
        f.write("# current      = CC\n")
        f.write("# units_sigma  = cm^2\n")
        f.write("# units_E      = GeV\n")
        f.write(f"# M_Z_GeV      = {TC.MZ}\n")
        f.write(f"# dipole_model = nenhum; colinear NLO, yadism {sf.yadism},"
                f" {sf.pdf_set}\n")
        f.write("# generated_by = tools/make_collinear_table.py (PHASIS)"
                f" commit {commit}\n")
        f.write(f"# M_W_GeV      = {TC.MW}\n")
        f.write(f"# pdf_set      = {sf.pdf_set}\n")
        f.write(f"# PTO          = {sf.pto} (0=LO, 1=NLO)\n")
        f.write(f"# esquema      = ZM-VFNS, TMC desligada, CKM completa\n")
        f.write(f"# formula      = mesma de dipole/src/sigma_nuN_core.cpp,"
                f" sem o fator (1-x)^7\n")
        f.write(f"# Q2min_GeV2   = {a.q2min}\n")
        f.write(f"# nos_por_decada_Q = {a.nodesQ}   nos_por_decada_x = {a.nodesX}\n")
        f.write(f"# tabela_SF    = {sf.lx.size} x {sf.lq.size} nos em (x, Q^2);"
                f" x em [{sf.x_min:.3g}, {sf.x_max:.3g}],"
                f" Q^2 em [{sf.q2_min:.3g}, {sf.q2_max:.3g}]\n")
        f.write("#\n")
        f.write("# AVISO -- EXTRAPOLACAO DO PDF. Abaixo de x = "
                f"{sf.pdf_xmin:g} o LHAPDF extrapola (continuation), e a\n")
        f.write("# ultima coluna diz que fracao de sigma vem dali. Nao e um\n")
        f.write("# defeito escondido: e exatamente a hipotese que o modelo de\n")
        f.write("# dipolo substitui. Onde essa fracao e grande, a inclinacao\n")
        f.write("# alpha_eff colinear e uma PROPRIEDADE DA EXTRAPOLACAO (uma\n")
        f.write("# lei de potencia continuada), nao uma medida.\n")
        f.write("# Enu_GeV sigma_GeV_minus2 sigma_cm2 fracao_x_abaixo_do_grid\n")
        for E, sv, fr in zip(Es, sig, frac):
            f.write(f"{E:.10g} {sv:.10g} {sv*TC.GEV2_TO_CM2:.10g} {fr:.6g}\n")
    print(f"Arquivo gerado: {a.out}")


if __name__ == "__main__":
    main()
