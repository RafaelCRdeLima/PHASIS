#ifndef PHASIS_DENSITY_HPP
#define PHASIS_DENSITY_HPP

#include <string>

namespace phasis {

// =====================================================================
// Perfil de materia.
//
// rho(r, theta) em g/cm^3. theta e o angulo polar usual, medido do eixo
// z; theta = pi/2 e o plano equatorial.
//
// O par (r_support_min, r_support_max) delimita onde rho pode ser nao
// nulo. Nao e cosmetico: e o que permite ao integrador cortar a faixa
// vazia em vez de pedir a quadratura adaptativa que descubra sozinha
// uma descontinuidade em degrau -- o que ela faz mal e caro.
// =====================================================================
struct DensityProfile {
    virtual ~DensityProfile() = default;

    virtual double rho(double r, double theta) const = 0;

    virtual double r_support_min() const { return 0.0; }
    virtual double r_support_max() const = 0;

    // Se true, rho nao depende de theta e o integrador pode pular o
    // calculo do angulo ao longo do raio (uma quadratura aninhada).
    // A Fase 3 introduz perfis com is_spherical() == false.
    virtual bool is_spherical() const { return true; }

    virtual std::string name() const = 0;
};

// ---------------------------------------------------------------------
// Esfera homogenea: rho = rho0 para r < R.
class UniformBall final : public DensityProfile {
public:
    UniformBall(double rho0, double R) : rho0_(rho0), R_(R) {}

    double rho(double r, double) const override {
        return (r < R_) ? rho0_ : 0.0;
    }
    double r_support_max() const override { return R_; }
    std::string name() const override { return "UniformBall"; }

    double rho0() const { return rho0_; }
    double R() const { return R_; }

private:
    double rho0_;
    double R_;
};

// ---------------------------------------------------------------------
// Halo em lei de potencia: rho = rho0 * (r/r0)^(-p),  r_in <= r <= r_out.
class PowerLawHalo final : public DensityProfile {
public:
    PowerLawHalo(double rho0, double r0, double p, double r_in, double r_out)
        : rho0_(rho0), r0_(r0), p_(p), r_in_(r_in), r_out_(r_out) {}

    double rho(double r, double) const override;

    double r_support_min() const override { return r_in_; }
    double r_support_max() const override { return r_out_; }
    std::string name() const override { return "PowerLawHalo"; }

    double rho0() const { return rho0_; }
    double r0()   const { return r0_; }
    double p()    const { return p_; }
    double r_in() const { return r_in_; }
    double r_out()const { return r_out_; }

private:
    double rho0_, r0_, p_, r_in_, r_out_;
};

} // namespace phasis

#endif
