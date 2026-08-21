#ifndef DIPOLE_MODELS_HPP
#define DIPOLE_MODELS_HPP

namespace dipole {

struct GBWParameters {
    double sigma0_mb = 29.12;
    double lambda = 0.277;
    double x0 = 0.41e-4;
    double Q0sq = 1.0;
};

struct IIMParameters {
    double BCGC = 5.5;
    double gamma_s = 0.6492;
    double N0 = 0.3658;
    double x0 = 0.00069e-6;
    double lambda = 0.2023;

    double bMax = 20.0;
    int Nb = 80;  //80 para desenvovlimento e 400 para versão final
};

double mbToGeVminus2(double sigma_mb);

double R0sq_GBW(double x, const GBWParameters& p);
double sigmaDipoleGBW(double r, double x, const GBWParameters& p);

double Qs_bCGC(double x, double b, const IIMParameters& p);
double amplitude_bCGC(double r, double x, double b, const IIMParameters& p);
double sigmaDipoleIIM(double r, double x, const IIMParameters& p);

} // namespace dipole




#endif