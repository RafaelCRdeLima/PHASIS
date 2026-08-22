#ifndef PHASIS_RAY_HPP
#define PHASIS_RAY_HPP

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
//   psi_node_rad     angulo, dentro do plano, entre o nodo ascendente
//                    e o ponto de retorno
//
// Um ponto do raio a angulo interno psi (medido do nodo) tem
//   cos(theta) = sin(inclination) * sin(psi)
//
// Com inclination = 0 isso da theta = pi/2 para todo psi: raio
// equatorial. E o default, e o unico caso que a Fase 1 exercita, ja que
// todos os perfis desta fase sao esfericos.
// =====================================================================
struct Ray {
    double E_inf_GeV       = 0.0;   // energia conservada, medida no infinito
    double b_cm            = 0.0;   // parametro de impacto b = L / E_inf
    double inclination_rad = 0.0;
    double psi_node_rad    = 0.0;
};

struct Result {
    double tau              = 0.0;
    double P_surv           = 1.0;
    double r_min_cm         = 0.0;   // ponto de retorno
    double column_density   = 0.0;   // X = integral de rho dl, em g/cm^2
    double path_length_cm   = 0.0;   // integral de dl dentro do suporte

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

    long   n_evals          = 0;
};

} // namespace phasis

#endif
