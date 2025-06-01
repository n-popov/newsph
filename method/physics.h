#pragma once

#include "../config/simulation_config.h"
#include "../utils/helpers.h"

void compute_eos_stiffened_gas(
    mysph::Particle<double>& particle,
    double gamma1, double C1, double ro1, 
    double gamma2, double C2, double ro2  
) {
    double gamma_eos, c0_eos, rho0_eos;

    if (particle.material == 0) { 
        gamma_eos = gamma1;
        c0_eos    = C1;
        rho0_eos  = ro1;
    } else { 
        gamma_eos = gamma2;
        c0_eos    = C2;
        rho0_eos  = ro2;
    }

    particle.p = c0_eos * c0_eos * (particle.rho - rho0_eos) + \
                 (gamma_eos - 1.0) * particle.rho * particle.e;

    if (particle.rho > 1e-9) { 
        double term_in_sqrt = c0_eos * c0_eos +
                              (gamma_eos - 1.0) * (particle.e + particle.p / particle.rho);
        if (term_in_sqrt >= 0) {
            particle.cs = std::sqrt(term_in_sqrt);
        } else {
            particle.cs = c0_eos;
        }
    } else {
        particle.cs = c0_eos; 
    }
}

void compute_velocity_gradient(std::vector<mysph::Particle<double>>& particles, int i,
                               const std::vector<mysph::Particle<double>>& neighbors,
                               double h) { 
    auto& pi = particles[i];
    pi.v00 = 0.0;
    pi.v01 = 0.0;
    pi.v02 = 0.0;
    pi.v10 = 0.0;
    pi.v11 = 0.0;
    pi.v12 = 0.0;
    pi.v20 = 0.0;
    pi.v21 = 0.0;
    pi.v22 = 0.0;

    for (const auto& pj : neighbors) {
        auto rij = pi.r - pj.r;
        auto vij_corrected = pj.v - pi.v; 
        auto grad_k = mysph::grad_kernel(rij, h); 

        double factor = pj.m / pj.rho;

        pi.v00 += factor * vij_corrected[0] * grad_k[0];
        pi.v01 += factor * vij_corrected[0] * grad_k[1];
        pi.v02 += factor * vij_corrected[0] * grad_k[2];

        pi.v10 += factor * vij_corrected[1] * grad_k[0];
        pi.v11 += factor * vij_corrected[1] * grad_k[1];
        pi.v12 += factor * vij_corrected[1] * grad_k[2];

        pi.v20 += factor * vij_corrected[2] * grad_k[0];
        pi.v21 += factor * vij_corrected[2] * grad_k[1];
        pi.v22 += factor * vij_corrected[2] * grad_k[2];
    }
}

// Von Mises plasticity model
void apply_von_mises_plasticity(mysph::Particle<double>& p) {
    // Calculate stress deviator second invariant J2
    p.J2 = 0.5 * (p.s00*p.s00 + p.s11*p.s11 + p.s22*p.s22 + 
                   2*(p.s01*p.s01 + p.s02*p.s02 + p.s12*p.s12));
    
    // Check yield criterion
    if (p.J2 > (p.Yo * p.Yo / 3.0)) {
        // Calculate yield factor
        double yield_factor = p.Yo / (std::sqrt(3.0 * p.J2));
        
        // Scale down deviatoric stress components
        p.s00 *= yield_factor;
        p.s01 *= yield_factor;
        p.s02 *= yield_factor;
        p.s11 *= yield_factor;
        p.s12 *= yield_factor;
        p.s22 *= yield_factor;
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
    S[0][0] = p.s00 - p.p; S[0][1] = p.s01;        S[0][2] = p.s02;
    S[1][0] = p.s01;        S[1][1] = p.s11 - p.p; S[1][2] = p.s12;
    S[2][0] = p.s02;        S[2][1] = p.s12;        S[2][2] = p.s22 - p.p;

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


// Modified compute_stress_rate
void compute_stress_rate_and_artificial_terms(mysph::Particle<double>& p, double dt, double artificial_stress_eps) {
    // Calculate deviatoric strain rate
    double div_v = p.v00 + p.v11 + p.v22;
    double eps00 = p.v00 - div_v/3.0;
    double eps01 = 0.5 * (p.v01 + p.v10);
    double eps02 = 0.5 * (p.v02 + p.v20);
    double eps11 = p.v11 - div_v/3.0;
    double eps12 = 0.5 * (p.v12 + p.v21);
    double eps22 = p.v22 - div_v/3.0;

    // Update deviatoric stress using Hooke's law (elastic trial stress)
    p.s00 += 2 * p.G * eps00 * dt;
    p.s01 += 2 * p.G * eps01 * dt;
    p.s02 += 2 * p.G * eps02 * dt;
    p.s11 += 2 * p.G * eps11 * dt;
    p.s12 += 2 * p.G * eps12 * dt;
    p.s22 += 2 * p.G * eps22 * dt;

    // Apply plasticity (modifies p.s00, p.s01, etc.)
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

            double s00a_dev = pi.s00; double s01a_dev = pi.s01; double s02a_dev = pi.s02;
            double s11a_dev = pi.s11; double s12a_dev = pi.s12;
            double s22a_dev = pi.s22;

            double s00b_dev = pj.s00; double s01b_dev = pj.s01; double s02b_dev = pj.s02;
            double s11b_dev = pj.s11; double s12b_dev = pj.s12;
            double s22b_dev = pj.s22;

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
