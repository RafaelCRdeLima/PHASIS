#!/usr/bin/env python3
"""Condicionamento da quadratura de sigma_nuN.

sigma(E) e uma funcao suave: sigma ~ E^0.36, entao uma perturbacao relativa dE/E
deve produzir dsigma/sigma ~ 0.36 * dE/E. A amplificacao |dsigma/sigma| / |dE/E|
deve valer ~0.4, nunca mais que ~5.

Na quadratura original ela vale 1e2 a 1e5: os nos de Simpson em ln(Q^2) cobrem
ln(1) ate ln(0.999*2*M_N*E), entao mexer em E desloca TODOS os nos, e como o pico
do integrando em Q^2 ~ M_W^2 e resolvido por apenas ~2-5 nos, o resultado pula.

E a mesma causa dos 29 passos decrescentes na tabela publicada. Este e o teste de
aceitacao da fase F4 e o diagnostico mais barato de saude da quadratura: nao
precisa de valor de referencia externo, so da suavidade que a fisica exige.

    python3 tests/test_conditioning.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "oracle"))
import dipole_oracle as o  # noqa: E402

# amplificacao maxima aceitavel; a fisica da ~0.4
MAX_AMPLIFICACAO = 5.0

ENERGIAS = (1.0e5, 1.0e6, 1.0e7, 1.0e9, 1.0e11, 1.0e13)
PERTURBACOES = (1.0e-7, 1.0e-6, 1.0e-5)


def main(**quad) -> int:
    print("Condicionamento da quadratura")
    print("  fisica: sigma ~ E^0.36  =>  amplificacao esperada ~0.4")
    print("  limite aceito: %.1f\n" % MAX_AMPLIFICACAO)
    print("  %-12s %-10s %-13s %s" % ("E [GeV]", "dE/E", "dsigma/sigma", "amplificacao"))

    pior = 0.0
    for E in ENERGIAS:
        s0 = o.sigma_cm2(E, **quad)
        for d in PERTURBACOES:
            s1 = o.sigma_cm2(E * (1.0 + d), **quad)
            rel = (s1 - s0) / s0
            amp = abs(rel / d)
            pior = max(pior, amp)
            flag = "" if amp <= MAX_AMPLIFICACAO else "   <-- FALHA"
            print("  %-12.2e %-10.0e %+-13.3e %-10.0f%s" % (E, d, rel, amp, flag))
        print()

    print("  pior amplificacao: %.4g  (limite %.1f)" % (pior, MAX_AMPLIFICACAO))
    if pior > MAX_AMPLIFICACAO:
        print("  RESULTADO: FALHA - sigma(E) nao e suave; a tabela carrega ruido de quadratura")
        return 1
    print("  RESULTADO: OK")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
