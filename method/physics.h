#pragma once

#include <stdexcept>

#include "../config/simulation_config.h"
#include "../utils/helpers.h"

void compute_eos_stiffened_gas(
    mysph::Particle<double>& particle,
    const config::SimulationConfig& config
) {
    auto props = get_material_properties(particle, config);

    const double& gamma = props.gruneisen_gamma;
    const double& c0 = props.sound_speed_coefficient;
    const double& rho0 = props.density;

    particle.p = c0 * c0 * (particle.rho - rho0) +
        (gamma - 1.0) * particle.rho * particle.e;

    double cs_squared = c0 * c0 +
                          (gamma - 1.0) * (particle.e + particle.p / particle.rho);

    particle.cs = std::sqrt(cs_squared);
}

void compute_velocity_gradient(mysph::Particle<double>& pi,
                               const std::vector<mysph::Particle<double>*>& neighbors,
                               double h) { 
    pi.v_grad = {};

    for (const auto ppj : neighbors) {
        const auto& pj = *ppj;
        auto rij = pi.r - pj.r;
        auto vij = pi.v - pj.v; 

        pi.v_grad = pi.v_grad - (pj.m / pj.rho) * matmul(vij, mysph::grad_kernel(rij, h));
    }
}

void apply_von_mises_plasticity(mysph::Particle<double>& p) {
    // Stress deviator second invariant J2
    p.J2 = 0.5 * (p.stress * p.stress);
    
    // Yield criterion
    if (p.J2 > (p.Yo * p.Yo / 3.0)) {
        double yield_factor = p.Yo / (std::sqrt(3.0 * p.J2));

        p.stress = p.stress * yield_factor;
    }
}

// Function to compute Monaghan Artificial Stress
void compute_monaghan_artificial_stress(mysph::Particle<double>& p, double eps) {
    double rhoi = p.rho;
    if (rhoi < 1e-9) {
        p.a_stress = {};

        return;
    }
    double rhoi21 = 1.0 / (rhoi * rhoi);

    mysph::vec3<double> rd;     // Artificial stress components in principal directions

    // Total stress tensor S = s_dev - p * I
    auto S = p.stress - mysph::SINGILAR9 * p.p;

    // Compute the principal stresses and eigenvectors
    auto [R_mat, V] = eigen_decomposition_3x3(S); // Fills R_mat - eigenvectors (columns) and V – eigenvalues, principal stresses

    // Artificial stress corrections in principal directions
    rd[0] = (V[0] > 0) ? -eps * V[0] * rhoi21 : 0.0;
    rd[1] = (V[1] > 0) ? -eps * V[1] * rhoi21 : 0.0;
    rd[2] = (V[2] > 0) ? -eps * V[2] * rhoi21 : 0.0;

    // Transform artificial stresses back to original frame: Rab = R_mat * diag(rd) * R_mat_transpose
    p.a_stress = transform_diag_inv_3x3(rd, R_mat); // Artificial stress tensor in original coordinates
}

void compute_stress_rate_and_artificial_terms(mysph::Particle<double>& p, double dt, double artificial_stress_eps) {
    // deviatoric strain rate
    auto strain_rate = 0.5 * (p.v_grad + transpose(p.v_grad)) - (div(p.v_grad) / 3.) * mysph::SINGILAR9;

    // Hooke's law (elastic trial stress)
    p.stress = (2 * p.G * dt) * strain_rate;

    apply_von_mises_plasticity(p);

    if (artificial_stress_eps > 0) { 
        compute_monaghan_artificial_stress(p, artificial_stress_eps);
    } else {
        p.a_stress = {};
    }
}

void compute_force(mysph::Particle<double>& pi, const std::vector<mysph::Particle<double>*>& neighbors, const config::SPHParameters& sph_params) {
    pi.F = {0.0, 0.0, 0.0};
            
    for (auto ppj : neighbors) {
        const auto& pj = *ppj;
        auto rij = pi.r - pj.r;
        auto vij = pi.v - pj.v;
        double r = mysph::abs(rij);

        if (r > 1e-9 * sph_params.h) {
            mysph::vector3d DWIJ = mysph::grad_kernel(rij, sph_params.h);
            double WIJ = mysph::kernel(rij, sph_params.h);

            double pa = pi.p;
            double pb = pj.p;
            double rhoa = pi.rho;
            double rhob = pj.rho;

            double rhoa21 = 1.0 / (rhoa * rhoa);
            double rhob21 = 1.0 / (rhob * rhob);

            // total stress tensor
            auto sa = pi.stress - mysph::SINGILAR9 * pa;
            auto sb = pj.stress - mysph::SINGILAR9 * pb;

            // artificial stress tensor
            auto r_ab = pi.a_stress + pj.a_stress;

            auto wdeltap = mysph::kernel(mysph::vec3<double>{r, r, r}, sph_params.h);
            auto n_art_stress = 2; // TODO move to config

            if (wdeltap > 1e-9) {
                r_ab = r_ab * std::pow(WIJ / wdeltap, n_art_stress);
            } else {
                r_ab = {};
            }

            pi.F = pi.F + pi.m * pj.m * matmul(DWIJ, sa * rhoa21 + sb * rhob21 + r_ab);
        }
    }
}
