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
                               const std::vector<mysph::Particle<double>>& neighbors,
                               double h) { 
    pi.v_grad = {};

    for (const auto& pj : neighbors) {
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
    if (rhoi < 1e-9) { // Avoid division by zero if density is tiny
        p.as00 = 0.0; p.as01 = 0.0; p.as02 = 0.0;
        p.as11 = 0.0; p.as12 = 0.0; p.as22 = 0.0;
        return;
    }
    double rhoi21 = 1.0 / (rhoi * rhoi);

    double S[3][3]; // Total stress tensor
    double rd[3];     // Artificial stress components in principal directions

    // Construct total stress tensor S = s_dev - p * I
    // Assuming p.pressure is positive for compression
    S[0][0] = p.stress[0] - p.p;  S[0][1] = p.stress[1];        S[0][2] = p.stress[2];
    S[1][0] = p.stress[3];        S[1][1] = p.stress[4] - p.p;  S[1][2] = p.stress[5];
    S[2][0] = p.stress[6];        S[2][1] = p.stress[7];        S[2][2] = p.stress[8] - p.p;

    // Compute the principal stresses and eigenvectors
    auto [R_mat, V] = eigen_decomposition_3x3(S); // Fills R_mat - eigenvectors (columns) and V – eigenvalues, principal stresses

    // Artificial stress corrections in principal directions
    rd[0] = (V[0] > 0) ? -eps * V[0] * rhoi21 : 0.0;
    rd[1] = (V[1] > 0) ? -eps * V[1] * rhoi21 : 0.0;
    rd[2] = (V[2] > 0) ? -eps * V[2] * rhoi21 : 0.0;

    // Transform artificial stresses back to original frame: Rab = R_mat * diag(rd) * R_mat_transpose
    auto Rab = transform_diag_inv_3x3(rd, R_mat); // Artificial stress tensor in original coordinates

    // Store the values in the particle (assuming p has r00, r01 etc. members)
    p.as00 = Rab[0][0];
    p.as01 = Rab[0][1]; // or Rab[1][0] due to symmetry
    p.as02 = Rab[0][2]; // or Rab[2][0]
    p.as11 = Rab[1][1];
    p.as12 = Rab[1][2]; // or Rab[2][1]
    p.as22 = Rab[2][2];
}

void compute_stress_rate_and_artificial_terms(mysph::Particle<double>& p, double dt, double artificial_stress_eps) {
    // Deviatoric strain rate
    auto strain_rate = 0.5 * (p.v_grad + transpose(p.v_grad)) - (div(p.v_grad) / 3.) * mysph::SINGILAR9;

    // Hooke's law (elastic trial stress)
    p.stress = (2 * p.G * dt) * strain_rate;

    apply_von_mises_plasticity(p);

    // Compute Monaghan Artificial Stress based on the updated deviatoric stresses and current pressure
    if (artificial_stress_eps > 0) { // Only compute if eps is non-zero
        compute_monaghan_artificial_stress(p, artificial_stress_eps);
    } else {
        // Ensure artificial stress terms are zero if not computed
        p.as00 = 0.0; p.as01 = 0.0; p.as02 = 0.0;
        p.as11 = 0.0; p.as12 = 0.0; p.as22 = 0.0;
    }
}

void compute_force(mysph::Particle<double>& pi, const std::vector<mysph::Particle<double>>& neighbors, const config::SPHParameters& sph_params) {
    pi.F = {0.0, 0.0, 0.0};
            
    for (auto& pj : neighbors) {
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

            double s00a_dev = pi.stress[0]; double s01a_dev = pi.stress[1]; double s02a_dev = pi.stress[2];
            double s11a_dev = pi.stress[4]; double s12a_dev = pi.stress[5];
            double s22a_dev = pi.stress[8];

            double s00b_dev = pj.stress[0]; double s01b_dev = pj.stress[1]; double s02b_dev = pj.stress[2];
            double s11b_dev = pj.stress[4]; double s12b_dev = pj.stress[5];
            double s22b_dev = pj.stress[8];

            double sigma00a = s00a_dev - pa; double sigma01a = s01a_dev;       double sigma02a = s02a_dev;
            double sigma10a = s01a_dev;       double sigma11a = s11a_dev - pa; double sigma12a = s12a_dev;
            double sigma20a = s02a_dev;       double sigma21a = s12a_dev;       double sigma22a = s22a_dev - pa;

            double sigma00b = s00b_dev - pb; double sigma01b = s01b_dev;       double sigma02b = s02b_dev;
            double sigma10b = s01b_dev;       double sigma11b = s11b_dev - pb; double sigma12b = s12b_dev;
            double sigma20b = s02b_dev;       double sigma21b = s12b_dev;       double sigma22b = s22b_dev - pb;

            double r00_ab = pi.as00 + pj.as00; double r01_ab = pi.as01 + pj.as01; double r02_ab = pi.as02 + pj.as02;
            double r11_ab = pi.as11 + pj.as11; double r12_ab = pi.as12 + pj.as12;
            double r22_ab = pi.as22 + pj.as22;

            double fab_pow_n = 0.0;
            double art_stress00 = 0.0, art_stress01 = 0.0, art_stress02 = 0.0;
            double art_stress11 = 0.0, art_stress12 = 0.0;
            double art_stress22 = 0.0;

            auto wdeltap = mysph::kernel(mysph::vector3d{r, r, r}, sph_params.h);
            auto n_art_stress = 2;

            if (wdeltap > 1e-9) {
                double fab = WIJ / wdeltap;
                fab_pow_n = std::pow(fab, n_art_stress);

                art_stress00 = fab_pow_n * r00_ab;
                art_stress01 = fab_pow_n * r01_ab;
                art_stress02 = fab_pow_n * r02_ab;
                art_stress11 = fab_pow_n * r11_ab;
                art_stress12 = fab_pow_n * r12_ab;
                art_stress22 = fab_pow_n * r22_ab;
            }

            double mb = pj.m;
            double ma = pi.m; 
            mysph::vector3d acc_contrib_vec;

            double term_xx = (sigma00a * rhoa21 + sigma00b * rhob21 + art_stress00);
            double term_xy = (sigma01a * rhoa21 + sigma01b * rhob21 + art_stress01);
            double term_xz = (sigma02a * rhoa21 + sigma02b * rhob21 + art_stress02);
            acc_contrib_vec[0] = mb * (term_xx * DWIJ[0] + term_xy * DWIJ[1] + term_xz * DWIJ[2]);

            double term_yx = (sigma10a * rhoa21 + sigma10b * rhob21 + art_stress01);
            double term_yy = (sigma11a * rhoa21 + sigma11b * rhob21 + art_stress11);
            double term_yz = (sigma12a * rhoa21 + sigma12b * rhob21 + art_stress12);
            acc_contrib_vec[1] = mb * (term_yx * DWIJ[0] + term_yy * DWIJ[1] + term_yz * DWIJ[2]);

            double term_zx = (sigma20a * rhoa21 + sigma20b * rhob21 + art_stress02);
            double term_zy = (sigma21a * rhoa21 + sigma21b * rhob21 + art_stress12);
            double term_zz = (sigma22a * rhoa21 + sigma22b * rhob21 + art_stress22);
            acc_contrib_vec[2] = mb * (term_zx * DWIJ[0] + term_zy * DWIJ[1] + term_zz * DWIJ[2]);

            pi.F = pi.F + ma * acc_contrib_vec;
        }
    }
}
