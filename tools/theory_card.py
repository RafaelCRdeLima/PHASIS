"""Cartao de teoria unico, compartilhado por todos os geradores colineares.

As constantes eletrofracas e a massa do nucleon sao COPIADAS de
PHASIS/dipole/include/parameters.hpp. Divergir aqui seria comparar maca
com laranja no lugar exato onde o teste T39 nao pegaria: uma diferenca
de normalizacao some em R(b,E), mas uma diferenca em M_W muda a FORMA,
que e justamente o que se pretende medir.
"""

# --- de dipole/include/parameters.hpp -------------------------------
GF      = 1.1663787e-5      # GeV^-2
MW      = 80.379            # GeV
MZ      = 91.1876           # GeV
MN      = 0.938272          # GeV
SIN2TW  = 0.23122
GEV2_TO_CM2 = 0.389379e-27

PDF_SET = "NNPDF31_nlo_as_0118"
PDF_MEMBER = 0
PDF_XMIN = 1.0e-9           # do .info; abaixo disso o LHAPDF EXTRAPOLA
PDF_QMIN = 1.65             # GeV
PDF_QMAX = 1.0e5            # GeV

PTO = 1                     # NLO


def theory(pto=PTO):
    """Cartao yadism/eko. NLO, ZM-VFNS, sem correcao de massa do alvo."""
    return dict(
        PTO=pto, PTODIS=pto, QED=0, FNS="ZM-VFNS", NfFF=5,
        ModEv="EXA", ModSV="expanded", XIF=1.0, XIR=1.0,
        Q0=1.65, nf0=4, Qref=MZ, nfref=5, alphas=0.118,
        alphaqed=0.007496252, MaxNfPdf=5, MaxNfAs=5,
        # limiares NNPDF3.1
        mc=1.51, mb=4.92, mt=172.5,
        kcThr=1.0, kbThr=1.0, ktThr=1.0,
        Qmc=1.51, Qmb=4.92, Qmt=172.5, HQ="POLE",
        MP=MN, MW=MW, MZ=MZ, GF=GF, SIN2TW=SIN2TW,
        CKM=("0.97428 0.22530 0.003470 "
             "0.22520 0.97345 0.041000 "
             "0.00862 0.04030 0.999152"),
        TMC=0, TMCM=0.0, n3lo_cf_variation=0, FONLLParts="full",
        IC=0, IB=0, DAMP=0, EScaleVar=1, ProtonMass=MN,
    )


def observable_card(points, xgrid):
    """CC, projetil nu (nao nubar), alvo isoescalar.

    Em yadism o observavel chamado 'F3' e xF3, nao F3 -- conferido em
    tools/check_yadism_convention.py contra a formula de ordem dominante.
    """
    return dict(
        interpolation_xgrid=list(xgrid),
        interpolation_is_log=True,
        interpolation_polynomial_degree=4,
        prDIS="CC",
        ProjectileDIS="neutrino",
        PolarizationDIS=0.0,
        PropagatorCorrection=0.0,
        TargetDIS="isoscalar",
        NCPositivityCharge=None,
        observables={
            "F2_total": points,
            "FL_total": points,
            "F3_total": points,
        },
    )
