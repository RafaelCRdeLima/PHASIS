#ifndef PHASIS_RAY_HPP
#define PHASIS_RAY_HPP

#include <string>

namespace phasis {

// =====================================================================
// Um raio nulo vindo do infinito.
//
// Em metrica esfericamente simetrica o movimento e PLANAR: o plano
// orbital contem a origem e e fixo. Basta entao girar o plano para
// cobrir qualquer geodesica -- nao ha necessidade de integrar em theta.
//
// A orientacao do plano entra por:
//   inclination_rad  angulo entre o plano orbital e o equatorial
// Eixo do disco = z. O plano orbital tem normal a angulo i do eixo z, e
// psi e a posicao azimutal DENTRO do plano orbital, medida da linha de
// nodos. A colatitude do ponto e entao
//
//     cos(theta) = sin(i) * sin(psi)
//
// Checagem: i = 0     -> theta = pi/2 sempre (raio no plano do disco);
//           i = pi/2  -> cos(theta) = sin(psi), cruza o disco em
//                        psi = 0 e psi = pi.
//
// psi_turn_rad e o valor de psi NO PONTO DE RETORNO. Ao longo do raio,
// psi cresce monotonicamente: o ramo de entrada vai de psi_in ate
// psi_turn, o de saida de psi_turn ate psi_out. Os dois ramos varrem
// intervalos de psi diferentes, e por isso deixam de ser iguais assim
// que rho depende de theta.
// =====================================================================
struct Ray {
    double E_inf_GeV       = 0.0;   // energia conservada, medida no infinito
    double b_cm            = 0.0;   // parametro de impacto b = L / E_inf
    double inclination_rad = 0.0;   // i
    double psi_turn_rad    = 0.0;   // psi no ponto de retorno

    // Emissao a distancia finita. r_emit_cm <= 0 significa "vindo do
    // infinito", que e o comportamento das Fases 1-4 e continua sendo o
    // default. Ver emission.hpp.
    double r_emit_cm = 0.0;
    bool   outward   = false;       // direcao no ponto de emissao
};

struct Result {
    double tau              = 0.0;
    double P_surv           = 1.0;
    double r_min_cm         = 0.0;   // ponto de retorno
    double column_density   = 0.0;   // X = integral de rho dl, em g/cm^2
    double path_length_cm   = 0.0;   // integral de dl dentro do suporte

    double tau_inbound      = 0.0;   // ramo de entrada, psi_turn - psi_off
    double tau_outbound     = 0.0;   // ramo de saida,   psi_turn + psi_off
    double column_inbound   = 0.0;   // g/cm^2
    double column_outbound  = 0.0;

    double theta_min_rad    = 0.0;   // extremos de colatitude visitados
    double theta_max_rad    = 0.0;
    int    n_disk_crossings = 0;     // vezes que o raio cruza z = 0

    double deflection_rad   = 0.0;   // Delta_phi = 2*INT(dphi/dr)dr - pi
    double E_loc_max_GeV    = 0.0;   // maior energia local ao longo do trecho util
    double winding_turns    = 0.0;   // (Delta_phi + pi) / 2pi

    bool   crosses_matter   = false;
    bool   captured         = false;
    bool   near_critical    = false; // T(r_t)/r_t <= 1e-6: raio enrola na
                                     // esfera de fotons; tau diverge
                                     // LOGARITMICAMENTE em b -> b_crit+.
                                     // Isso e fisico, nao falha numerica.
    bool   tolerance_met    = true;  // false se a quadratura esgotou profundidade

    bool   error            = false; // excecao capturada na varredura
    std::string error_msg;

    long   n_evals          = 0;
};

} // namespace phasis

#endif
