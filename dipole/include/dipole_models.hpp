#ifndef DIPOLE_MODELS_HPP
#define DIPOLE_MODELS_HPP

namespace dipole {

// Golec-Biernat & Wusthoff, Phys. Rev. D59 (1999) 014017 / D60 (1999) 114023.
// Ajuste de QUATRO sabores (inclui charm), o mesmo usado por
// Kutak & Kwiecinski EPJ C29 (2003) 521.
//
// As massas efetivas fazem parte do ajuste -- nao sao convencao livre.
struct GBWParameters {
    double sigma0_mb = 29.12;
    double lambda = 0.277;
    double x0 = 0.41e-4;
    double Q0sq = 1.0;

    double m_light = 0.14;   // u, d, s
    double m_charm = 1.5;
};

// b-CGC.  Rezaeian & Schmidt, Phys. Rev. D 88 (2013) 074016, Tabela II,
// segunda linha (m_c = 1.4 GeV):
//
//   B_CGC = 5.5 GeV^-2   gamma_s = 0.6492 +- 0.0003   N0 = 0.3658 +- 0.0006
//   x0 = 0.00069 +- 6.46e-6   lambda = 0.2023 +- 0.0003   chi2/dof = 1.249
//
// ATENCAO ao x0: o valor central e 0.00069 = 6.9e-4, e 6.46e-6 e a
// INCERTEZA. Este codigo trazia x0 = 0.00069e-6 = 6.9e-10, ou seja um
// fator 1e6 a menos, por ter lido o expoente da incerteza como
// multiplicador do valor central. Com isso o modelo dava <mod/dado> =
// 0.10 contra o HERA -- o dado a que ele foi ajustado. Ver
// CAMPANHA_CORRECAO.md secao 4, Q2.
//
// O artigo usa massas leves m_u = 1e-2 a 1e-4 GeV (praticamente sem
// massa) e m_c = 1.4 GeV; kappa = 9.9 e o valor LO do BFKL; e a
// rapidez e Y = ln(1/x).
struct IIMParameters {
    double BCGC = 5.5;
    double gamma_s = 0.6492;
    double N0 = 0.3658;
    double x0 = 0.69e-3;
    double lambda = 0.2023;

    double m_light = 0.01;   // artigo: 1e-2 a 1e-4 GeV
    double m_charm = 1.4;

    double bMax = 20.0;
    int Nb = 80;   // converge com 40; ver CAMPANHA_CORRECAO.md Q2
};

double mbToGeVminus2(double sigma_mb);

double R0sq_GBW(double x, const GBWParameters& p);
double sigmaDipoleGBW(double r, double x, const GBWParameters& p);

double Qs_bCGC(double x, double b, const IIMParameters& p);
double amplitude_bCGC(double r, double x, double b, const IIMParameters& p);
double sigmaDipoleIIM(double r, double x, const IIMParameters& p);

} // namespace dipole




#endif