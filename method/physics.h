#pragma once

#include <stdexcept>

#include "../config/simulation_config.h"
#include "../utils/helpers.h"

// verified
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

// verified
void compute_velocity_gradient(mysph::Particle<double>& pi,
                               const std::vector<mysph::Particle<double>*>& neighbors,
                               const config::SPHParameters& sph_params) { 
    pi.v_grad = {};

    for (const auto ppj : neighbors) {
        const auto& pj = *ppj;
        auto rij = pi.r - pj.r;
        auto vji = pj.v - pi.v; 

        pi.v_grad = pi.v_grad + (pj.m / pj.rho) * matmul(vji, mysph::grad_kernel(rij, sph_params.h, sph_params.kernel));
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

// validated
// Function to compute Monaghan Artificial Stress
void compute_monaghan_artificial_stress(mysph::Particle<double>& p, double eps) {
    if (eps <= 0.) {
        p.a_stress = {};

        return;
    }

    double rhoi = p.rho;
    double rhoi21 = 1.0 / (rhoi * rhoi);

    mysph::vec9<double> rd = {};     // Artificial stress in principal directions

    auto S = p.stress - mysph::SINGILAR9 * p.p;

    auto [R_mat, V] = eigen_decomposition_3x3(S);

    rd[0] = (V[0] > 0) ? -eps * V[0] * rhoi21 : 0.0;
    rd[4] = (V[1] > 0) ? -eps * V[1] * rhoi21 : 0.0;
    rd[8] = (V[2] > 0) ? -eps * V[2] * rhoi21 : 0.0;

    p.a_stress = matmul(matmul(R_mat, rd), transpose(R_mat));
}


void continuity_equation(mysph::Particle<double>& p, const std::vector<mysph::Particle<double>*>& neighbors, const config::SPHParameters& sph_params) {
    p.arho = 0.;

    for(auto ppb: neighbors) {
        auto& pb = *ppb;

        p.arho += pb.m + (p.v - pb.v) * mysph::grad_kernel(p.r - pb.r, sph_params.h, Kernel::TUTORIAL);
    }
}


// validated
void compute_hookes_deviatoric_stress_rate(mysph::Particle<double>& p) {
    auto epsilon = 0.5 * (p.v_grad + transpose(p.v_grad)) - (div(p.v_grad) / 3.) * mysph::SINGILAR9;
    auto omega = 0.5 * (p.v_grad - transpose(p.v_grad));

    p.acc_stress = 2 * p.G * epsilon + matmul(p.stress, transpose(omega)) + matmul(omega, transpose(p.stress));
}

// validated
void compute_force(mysph::Particle<double>& pi, const std::vector<mysph::Particle<double>*>& neighbors, const config::SPHParameters& sph_params) {
    pi.a = {0.0, 0.0, 0.0};
            
    for (auto ppj : neighbors) {
        const auto& pj = *ppj;
        auto rij = pi.r - pj.r;
        auto vij = pi.v - pj.v;
        double r = mysph::abs(rij);

        if (r > 1e-9 * sph_params.h) {
            mysph::vector3d DWIJ = mysph::grad_kernel(rij, sph_params.h, sph_params.kernel);
            double WIJ = mysph::kernel(rij, sph_params.h, sph_params.kernel);
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
            auto wdeltap = mysph::kernel(mysph::vec3<double>{r, r, r}, sph_params.h, sph_params.kernel);

            if (wdeltap > 1e-9) {
                r_ab = r_ab * std::pow(WIJ / wdeltap, sph_params.n_art_stress);
            } else {
                r_ab = {};
            }
            pi.a = pi.a + pj.m * matmul(DWIJ, sa * rhoa21 + sb * rhob21 + r_ab);
        }
    }
}

// validated
void compute_energy_with_stress(mysph::Particle<double>& pi, const std::vector<mysph::Particle<double>*>& neighbors, const config::SPHParameters& sph_params) {
    pi.ae = 0.;
    auto& h = sph_params.h;

    for (auto ppj: neighbors) {
        const auto& pj = *ppj;

        auto vij = pi.v - pj.v;
        auto rij = pi.r - pj.r;
        auto cij = (pi.cs + pj.cs) / 2;
        auto rhoij = (pi.rho + pj.rho) / 2;

        auto Piij = 0.;

        if (vij * rij < 0.) {
            auto phiij = (vij * rij * h) / (rij * rij + sph_params.avisc_eta * sph_params.avisc_eta * h * h);

            Piij = -sph_params.avisc_alpha * cij * phiij + sph_params.avisc_beta * phiij * phiij;
            Piij *= rhoij;
        }
        

        pi.ae += 0.5 * pj.m * (pi.p / (pi.rho * pi.rho) + pj.p / (pj.rho * pj.rho) + Piij); //* (vij * mysph::grad_kernel(rij, h, sph_params.kernel));
    }

    auto strain_rate = 0.5 * (pi.v_grad + transpose(pi.v_grad));
    pi.ae += (pi.stress * strain_rate) / pi.rho;
}