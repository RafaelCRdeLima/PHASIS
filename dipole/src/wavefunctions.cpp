#include "wavefunctions.hpp"

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>

namespace dipole {

double epsilon2(double z, double Q2, const Parameters& p)
{
    return z*(1.0 - z)*Q2
         + z*p.m*p.m
         + (1.0 - z)*p.mu*p.mu;
}

double psiT2_pre(double z, double Q2, const Parameters& p,
                 double eps2v, double K0, double K1)
{
    (void)Q2;

    const double gV2 = p.gV*p.gV;
    const double gA2 = p.gA*p.gA;

    const double zbar = 1.0 - z;

    const double termK1 =
        (gV2 + gA2)
        * (z*z + zbar*zbar)
        * eps2v
        * K1*K1;

    const double termK0 =
        (
            gV2 * std::pow(z*p.m + zbar*p.mu, 2)
          + gA2 * std::pow(z*p.m - zbar*p.mu, 2)
        )
        * K0*K0;

    const double prefactor =
        6.0 * p.alphaEW / std::pow(2.0*pi, 2);

    return prefactor * (termK1 + termK0);
}

double psiL2_pre(double z, double Q2, const Parameters& p,
                 double eps2v, double K0, double K1)
{
    const double gV2 = p.gV*p.gV;
    const double gA2 = p.gA*p.gA;

    const double zbar = 1.0 - z;

    const double m  = p.m;
    const double mu = p.mu;

    const double A =
        gV2 * std::pow(m - mu, 2)
      + gA2 * std::pow(m + mu, 2);

    const double Bv =
        2.0*Q2*z*zbar
      + (m - mu)*(z*m - zbar*mu);

    const double Ba =
        2.0*Q2*z*zbar
      + (m + mu)*(z*m + zbar*mu);

    const double termK1 =
        A * eps2v * K1*K1;

    const double termK0 =
        (
            gV2 * Bv*Bv
          + gA2 * Ba*Ba
        )
        * K0*K0;

    const double prefactor =
        6.0 * p.alphaEW / (std::pow(2.0*pi, 2) * Q2);

    return prefactor * (termK1 + termK0);
}

double psiT2(double r, double z, double Q2, const Parameters& p)
{
    const double eps2v = epsilon2(z, Q2, p);
    const double arg   = std::sqrt(eps2v)*r;

    return psiT2_pre(z, Q2, p, eps2v,
                     boost::math::cyl_bessel_k(0, arg),
                     boost::math::cyl_bessel_k(1, arg));
}

double psiL2(double r, double z, double Q2, const Parameters& p)
{
    const double eps2v = epsilon2(z, Q2, p);
    const double arg   = std::sqrt(eps2v)*r;

    return psiL2_pre(z, Q2, p, eps2v,
                     boost::math::cyl_bessel_k(0, arg),
                     boost::math::cyl_bessel_k(1, arg));
}

} // namespace dipole
