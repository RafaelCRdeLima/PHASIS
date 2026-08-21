#include "integrals.hpp"

#include "wavefunctions.hpp"

#include <cmath>

namespace dipole {

double simpson(
    const std::function<double(double)>& f,
    double a,
    double b,
    int N
)
{
    if (N % 2 != 0) {
        ++N;
    }

    const double h = (b - a)/N;

    double sum = f(a) + f(b);

    for (int i = 1; i < N; ++i) {

        const double x = a + i*h;

        if (i % 2 == 0) {
            sum += 2.0*f(x);
        } else {
            sum += 4.0*f(x);
        }
    }

    return sum*h/3.0;
}

double integrate2D(
    const std::function<double(double,double)>& f,
    double ax,
    double bx,
    int Nx,
    double ay,
    double by,
    int Ny
)
{
    auto gx = [&](double x)
    {
        auto gy = [&](double y)
        {
            return f(x,y);
        };

        return simpson(gy, ay, by, Ny);
    };

    return simpson(gx, ax, bx, Nx);
}

double integrandT(
    double r,
    double z,
    double x,
    double Q2,
    const Parameters& wf,
    const GBWParameters& gbw
)
{
    const double psi2 =
        psiT2(r, z, Q2, wf);

    const double sigma =
        sigmaDipoleGBW(r, x, gbw);

    return 2.0*pi*r*psi2*sigma;
}

double integrandL(
    double r,
    double z,
    double x,
    double Q2,
    const Parameters& wf,
    const GBWParameters& gbw
)
{
    const double psi2 =
        psiL2(r, z, Q2, wf);

    const double sigma =
        sigmaDipoleGBW(r, x, gbw);

    return 2.0*pi*r*psi2*sigma;
}

double FT_GBW(
    double x,
    double Q2,
    const Parameters& wf,
    const GBWParameters& gbw,
    int Nr,
    int Nz
)
{
    auto f = [&](double r, double z)
    {
        return integrandT(
            r, z, x, Q2, wf, gbw
        );
    };

    const double integral =
        integrate2D(
            f,
            wf.rMin, wf.rMax, Nr,
            wf.zMin, wf.zMax, Nz
        );

    return Q2/(4.0*pi*pi)*integral;
}

double FL_GBW(
    double x,
    double Q2,
    const Parameters& wf,
    const GBWParameters& gbw,
    int Nr,
    int Nz
)
{
    auto f = [&](double r, double z)
    {
        return integrandL(
            r, z, x, Q2, wf, gbw
        );
    };

    const double integral =
        integrate2D(
            f,
            wf.rMin, wf.rMax, Nr,
            wf.zMin, wf.zMax, Nz
        );

    return Q2/(4.0*pi*pi)*integral;
}

double F2_GBW(
    double x,
    double Q2,
    const Parameters& wf,
    const GBWParameters& gbw,
    int Nr,
    int Nz
)
{
    return
        FT_GBW(x,Q2,wf,gbw,Nr,Nz)
      + FL_GBW(x,Q2,wf,gbw,Nr,Nz);
}

double integrandT_IIM(
    double r,
    double z,
    double x,
    double Q2,
    const Parameters& wf,
    const IIMParameters& iim
)
{
    const double psi2 = psiT2(r, z, Q2, wf);
    const double sigma = sigmaDipoleIIM(r, x, iim);

    return 2.0*pi*r*psi2*sigma;
}

double integrandL_IIM(
    double r,
    double z,
    double x,
    double Q2,
    const Parameters& wf,
    const IIMParameters& iim
)
{
    const double psi2 = psiL2(r, z, Q2, wf);
    const double sigma = sigmaDipoleIIM(r, x, iim);

    return 2.0*pi*r*psi2*sigma;
}

double FT_IIM(
    double x,
    double Q2,
    const Parameters& wf,
    const IIMParameters& iim,
    int Nr,
    int Nz
)
{
    auto f = [&](double r, double z)
    {
        return integrandT_IIM(r, z, x, Q2, wf, iim);
    };

    const double integral =
        integrate2D(
            f,
            wf.rMin, wf.rMax, Nr,
            wf.zMin, wf.zMax, Nz
        );

    return Q2/(4.0*pi*pi)*integral;
}

double FL_IIM(
    double x,
    double Q2,
    const Parameters& wf,
    const IIMParameters& iim,
    int Nr,
    int Nz
)
{
    auto f = [&](double r, double z)
    {
        return integrandL_IIM(r, z, x, Q2, wf, iim);
    };

    const double integral =
        integrate2D(
            f,
            wf.rMin, wf.rMax, Nr,
            wf.zMin, wf.zMax, Nz
        );

    return Q2/(4.0*pi*pi)*integral;
}

double F2_IIM(
    double x,
    double Q2,
    const Parameters& wf,
    const IIMParameters& iim,
    int Nr,
    int Nz
)
{
    return
        FT_IIM(x,Q2,wf,iim,Nr,Nz)
      + FL_IIM(x,Q2,wf,iim,Nr,Nz);
}

} // namespace dipole