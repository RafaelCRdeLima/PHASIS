#ifndef PHASIS_DENSITY_HPP
#define PHASIS_DENSITY_HPP

#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

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
//
// INVARIANTE OBRIGATORIO -- rho tem de ser INCLUSIVA nas duas bordas:
//
//     rho(r_support_min) e rho(r_support_max) devem valer o LIMITE PELO
//     INTERIOR, nunca zero.
//
// Nao e detalhe de estilo. Os integradores usam exatamente [r_lo, r_hi]
// como intervalo, entao as bordas SAO pontos de avaliacao. Se rho for
// zero na borda e nao-zero um fio para dentro, o estagio de Runge-Kutta
// que cai ali ve um degrau, o estimador de erro nao converge e o passo
// colapsa ate h_min -- e encolher nunca ajuda, porque a descontinuidade
// esta no extremo.
//
// Custou dois bugs, um em cada sentido: na Fase 3 um trecho VAZIO
// tocava r_in (onde rho e corretamente nao-nula), resolvido nao
// integrando materia onde nao ha; na Fase 4 o ramo de entrada COMECAVA
// em r_out, onde UniformBall devolvia zero por usar `r < R` em vez de
// `r <= R`. Sao a mesma classe: borda dura num extremo de integracao.
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
        // `<=`, nao `<`: ver o invariante de borda em DensityProfile.
        return (r <= R_) ? rho0_ : 0.0;
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

// ---------------------------------------------------------------------
// Disco fino com abertura (flared), em coordenadas cilindricas
//
//     R = r sin(theta),   z = r cos(theta)
//
//     rho(R,z) = rho_0 (R/R_0)^(-p) exp( -z^2 / (2 H(R)^2) )
//     H(R)     = H_0 (R/R_0)^q
//
// zero fora de [r_in, r_out] no raio ESFERICO r.
//
// Defaults p = 15/8, q = 9/8: regiao externa de Shakura-Sunyaev.
//
// E o primeiro perfil com is_spherical() == false. A partir dele os dois
// ramos do raio deixam de ser iguais.
class FlaredThinDisk final : public DensityProfile {
public:
    FlaredThinDisk(double rho0, double R0, double H0,
                   double r_in, double r_out,
                   double p = 15.0/8.0, double q = 9.0/8.0)
        : rho0_(rho0), R0_(R0), H0_(H0), r_in_(r_in), r_out_(r_out), p_(p), q_(q) {}

    double rho(double r, double theta) const override;

    double r_support_min() const override { return r_in_; }
    double r_support_max() const override { return r_out_; }
    bool   is_spherical()  const override { return false; }
    std::string name() const override { return "FlaredThinDisk"; }

    double H_of_R(double R) const;

    double rho0() const { return rho0_; }
    double R0()   const { return R0_; }
    double H0()   const { return H0_; }
    double p()    const { return p_; }
    double q()    const { return q_; }

private:
    double rho0_, R0_, H0_, r_in_, r_out_, p_, q_;
};

// ---------------------------------------------------------------------
// ADAF quase-esferico: rho = rho_0 (r/R_0)^(-3/2), sem estrutura vertical.
//
// Ponte entre o caso esferico e o disco: mesma lei radial de um fluxo de
// acrecao dominado por adveccao, mas ainda com is_spherical() == true,
// entao serve de controle contra os perfis das Fases 1-2.
class QuasiSphericalADAF final : public DensityProfile {
public:
    QuasiSphericalADAF(double rho0, double R0, double r_in, double r_out)
        : rho0_(rho0), R0_(R0), r_in_(r_in), r_out_(r_out) {}

    double rho(double r, double) const override;

    double r_support_min() const override { return r_in_; }
    double r_support_max() const override { return r_out_; }
    std::string name() const override { return "QuasiSphericalADAF"; }

private:
    double rho0_, R0_, r_in_, r_out_;
};

// =====================================================================
// Registro de perfis, para o teste generico de borda.
//
// Duas ocorrencias do bug de borda em fases diferentes e assinatura de
// problema ESTRUTURAL, e comentario nao impede a terceira. O registro
// e AUTO-REGISTRAVEL: cada perfil se inscreve por um objeto estatico, e
// o teste varre o registro. Um perfil novo entra no teste sozinho, sem
// que ninguem lembre de adiciona-lo a uma lista.
// =====================================================================
class ProfileRegistry {
public:
    using Factory = std::function<std::shared_ptr<DensityProfile>()>;

    static ProfileRegistry& instance();

    void add(std::string nome, Factory f);
    const std::vector<std::pair<std::string, Factory>>& all() const { return itens_; }

private:
    std::vector<std::pair<std::string, Factory>> itens_;
};

struct ProfileRegistrar {
    ProfileRegistrar(std::string nome, ProfileRegistry::Factory f) {
        ProfileRegistry::instance().add(std::move(nome), std::move(f));
    }
};

#define PHASIS_REGISTER_PROFILE(tag, expr)                                  \
    static const ::phasis::ProfileRegistrar phasis_reg_##tag(              \
        #tag, []() -> std::shared_ptr<::phasis::DensityProfile> { return expr; })

} // namespace phasis

#endif
