"""Tabela colinear de F2, FL, xF3 para nu N -> l- X (CC, alvo isoescalar).

Roda yadism sobre uma grade (x, Q^2) e grava .npz.

MEMORIA -- a razao de existir o fatiamento
------------------------------------------
yadism.Runner.get_result() guarda a grade de funcoes de coeficiente de
CADA ponto pedido ate o fim da chamada. Medido nesta maquina:

    100 pontos -> 361 MB      600 pontos -> 837 MB
    ou seja  RSS ~ 266 MB + 0.95 MB por ponto

Pedir a grade inteira de uma vez (121 x 71 = 8591 pontos) custa ~8.4 GB.
Numa maquina de 15 GB com 2 GB de swap isso trava o sistema inteiro: nao
dispara o OOM-killer, so paginacao ate o desktop parar de responder.
Aconteceu.

Por isso o calculo vai FATIA POR FATIA em Q^2 (uma fatia = nx pontos,
~380 MB de pico) e cada fatia e gravada assim que sai. O consumo fica
CONSTANTE em relacao ao tamanho da grade, e uma interrupcao custa uma
fatia, nao a corrida inteira.

EXTRAPOLACAO -- registrada, nao escondida
----------------------------------------
Abaixo de x = 1e-9 e acima de Q = 1e5 GeV o LHAPDF extrapola. Nao e
defeito: e exatamente a hipotese que o dipolo substitui. A fracao de
sigma que vem dessa regiao e reportada por energia em
tools/make_collinear_table.py.
"""
import argparse
import logging
import os
import resource
import time
import warnings

import numpy as np

import theory_card as TC

logging.disable(logging.CRITICAL)
warnings.filterwarnings("ignore")


def rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nx", type=int, default=121)
    ap.add_argument("--nq", type=int, default=71)
    ap.add_argument("--xmin", type=float, default=1.0e-15)
    ap.add_argument("--xmax", type=float, default=1.0 - 1.0e-6)
    ap.add_argument("--q2min", type=float, default=1.0)
    ap.add_argument("--q2max", type=float, default=2.0e14)
    ap.add_argument("--pto", type=int, default=TC.PTO)
    ap.add_argument("--out", default="data/collinear_sf.npz")
    ap.add_argument("--work", default="data/collinear_sf_parcial")
    a = ap.parse_args()

    import lhapdf
    import yadism

    x = np.geomspace(a.xmin, a.xmax, a.nx)
    q2 = np.geomspace(a.q2min, a.q2max, a.nq)
    os.makedirs(a.work, exist_ok=True)

    # A base de interpolacao vai ABAIXO do menor x pedido: no no de borda
    # a propria base extrapolaria, e esse erro se somaria ao da
    # extrapolacao do PDF sem ser distinguivel dele.
    xgrid = np.geomspace(a.xmin*1.0e-1, 1.0, 120)
    theory = TC.theory(a.pto)
    pdf = lhapdf.mkPDF(TC.PDF_SET, TC.PDF_MEMBER)

    print(f"grade {a.nx} x {a.nq}; fatias de {a.nx} pontos em Q^2", flush=True)
    t0 = time.time()

    for iq, Q2 in enumerate(q2):
        fatia = os.path.join(a.work, f"q{iq:04d}.npz")
        if os.path.exists(fatia):
            continue   # retomada: fatia ja feita

        pts = [{"x": float(xx), "Q2": float(Q2), "y": 0.5} for xx in x]
        out = yadism.Runner(theory, TC.observable_card(pts, xgrid)).get_result()
        res = out.apply_pdf(pdf)
        del out   # solta a grade de coeficientes ANTES da proxima fatia

        np.savez(
            fatia,
            F2=np.array([p["result"] for p in res["F2_total"]]),
            FL=np.array([p["result"] for p in res["FL_total"]]),
            xF3=np.array([p["result"] for p in res["F3_total"]]),
        )
        del res

        dt = time.time() - t0
        feitas = iq + 1
        print(f"  Q^2 = {Q2:11.4g}  ({feitas}/{a.nq})  "
              f"RSS pico {rss_mb():.0f} MB  {dt:.0f} s  "
              f"restam ~{dt/feitas*(a.nq-feitas)/60:.0f} min", flush=True)

    F2 = np.empty((a.nq, a.nx))
    FL = np.empty((a.nq, a.nx))
    xF3 = np.empty((a.nq, a.nx))
    for iq in range(a.nq):
        d = np.load(os.path.join(a.work, f"q{iq:04d}.npz"))
        F2[iq], FL[iq], xF3[iq] = d["F2"], d["FL"], d["xF3"]

    np.savez_compressed(
        a.out, x=x, q2=q2, F2=F2, FL=FL, xF3=xF3,
        pto=a.pto, pdf_set=TC.PDF_SET, pdf_member=TC.PDF_MEMBER,
        pdf_xmin=TC.PDF_XMIN, pdf_qmax=TC.PDF_QMAX,
        yadism_version=yadism.__version__,
    )
    print(f"gravado {a.out}   RSS pico {rss_mb():.0f} MB", flush=True)
    print(f"  F2 negativo em {int((F2 < 0).sum())} nos de {F2.size}", flush=True)
    print(f"  FL negativo em {int((FL < 0).sum())} nos", flush=True)


if __name__ == "__main__":
    main()
