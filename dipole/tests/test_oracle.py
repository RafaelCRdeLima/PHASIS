#!/usr/bin/env python3
"""F0 - rede de seguranca: o oraculo Python reproduz o binario C++?

Roda sem LHAPDF e sem compilar nada. Compara contra o baseline congelado
em tests/baseline/sigma_baseline.json (gerado com --use-F3 0).

    python3 tests/test_oracle.py

Falha se o oraculo divergir do baseline em mais de 1e-4 relativo.
Se este teste falhar DEPOIS de uma mudanca no C++, o oraculo precisa ser
atualizado junto - deliberadamente, nunca por acidente.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "oracle"))
import dipole_oracle as o  # noqa: E402

TOL = 1.0e-4
ROOT = Path(__file__).resolve().parent.parent


def main() -> int:
    baseline = json.loads((ROOT / "tests/baseline/sigma_baseline.json").read_text())
    ref = baseline["GBW_semF3"]

    print("F0 - oraculo Python vs binario C++ (GBW, sem F3)")
    print("  %-14s %-14s %-14s %s" % ("E [GeV]", "oraculo", "baseline", "desvio"))

    worst = 0.0
    for log_e, expected in zip(ref["logE"], ref["sigma_cm2"]):
        E = 10.0 ** log_e
        got = o.sigma_cm2(E)
        rel = abs(got - expected) / expected
        worst = max(worst, rel)
        flag = "" if rel <= TOL else "  <-- FALHA"
        print("  %-14.6e %-14.6e %-14.6e %+.2e%s" % (E, got, expected, (got - expected) / expected, flag))

    print("\n  pior desvio relativo: %.2e  (tolerancia %.0e)" % (worst, TOL))
    if worst > TOL:
        print("  RESULTADO: FALHA")
        return 1
    print("  RESULTADO: OK - oraculo fiel, seguro para validar mudancas no C++")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
