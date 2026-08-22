#!/usr/bin/env python3
"""F0 - rede de seguranca: o oraculo Python continua fiel ao C++?

Duas verificacoes:

  1. MODO PRE-F2: com as grades uniformes originais, o oraculo reproduz
     tests/baseline/sigma_baseline.json - o sigma que o binario dava antes
     da campanha. Prova que o oraculo e uma reimplementacao fiel, e nao
     uma aproximacao que foi acompanhando as mudancas.

  2. MODO POS-F2: com a quadratura nova (log r + z nas duas pontas), o
     oraculo reproduz o valor de referencia de F_2 que o teste C++
     (make test-quadrature) mede.

Roda sem LHAPDF e sem compilar nada:

    python3 tests/test_oracle.py
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "oracle"))
import dipole_oracle as o  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent

TOL_BASELINE = 1.0e-4
TOL_F2 = 0.02

# grades uniformes originais, congeladas junto com o baseline
PRE_F2 = dict(Nr=50, Nz=30, logr=False, logz=False,
              rmin=1e-6, rmax=1e2, zmin=1e-6)

# valor de referencia de F_2(u,d) em x=1e-5, Q2=M_W^2 (ver src/test_quadrature.cpp)
F2_REF = 18.044
FL_F2_REF = 0.0986


def check_baseline() -> int:
    baseline = json.loads((ROOT / "tests/baseline/sigma_baseline.json").read_text())
    ref = baseline["GBW_semF3"]

    print("1) MODO PRE-F2 - oraculo vs binario congelado (GBW, sem F3)")
    print("   %-14s %-14s %-14s %s" % ("E [GeV]", "oraculo", "baseline", "desvio"))

    worst = 0.0
    for log_e, expected in zip(ref["logE"], ref["sigma_cm2"]):
        E = 10.0 ** log_e
        got = o.sigma_cm2(E, **PRE_F2)
        rel = abs(got - expected) / expected
        worst = max(worst, rel)
        flag = "" if rel <= TOL_BASELINE else "  <-- FALHA"
        print("   %-14.6e %-14.6e %-14.6e %+.2e%s"
              % (E, got, expected, (got - expected) / expected, flag))

    print("\n   pior desvio: %.2e  (tolerancia %.0e)  %s\n"
          % (worst, TOL_BASELINE, "OK" if worst <= TOL_BASELINE else "FALHA"))
    return 0 if worst <= TOL_BASELINE else 1


def check_pos_f2() -> int:
    print("2) MODO POS-F2 - F_2(u,d) em x=1e-5, Q2=M_W^2")

    FT, FL = o.F_TL(1.0e-5, o.MW ** 2, 0.03, 0.03, 200, 200)
    f2 = FT + FL
    ratio = FL / f2

    d2 = abs(f2 - F2_REF) / F2_REF
    dr = abs(ratio - FL_F2_REF) / FL_F2_REF

    print("   F_2     = %10.5f   (C++ %.3f, desvio %.2f%%)" % (f2, F2_REF, 100 * d2))
    print("   F_L/F_2 = %10.5f   (C++ %.4f, desvio %.2f%%)" % (ratio, FL_F2_REF, 100 * dr))

    ok = d2 <= TOL_F2 and dr <= TOL_F2
    print("\n   %s\n" % ("OK" if ok else "FALHA"))
    return 0 if ok else 1


def main() -> int:
    print("F0 - fidelidade do oraculo Python\n")
    falhas = check_baseline() + check_pos_f2()
    print("RESULTADO: %s" % ("OK" if falhas == 0 else "FALHA"))
    return 1 if falhas else 0


if __name__ == "__main__":
    raise SystemExit(main())
