#pragma once

#include "common.h"
#include "particle.h"

mysph::vec3<double> compute_xsph_corrected_velocities(
    const Particle<double>& pa,                     // All particles (destination and source)
    const std::vector<Particle<double>>& neighbors, // neighbors[i] is a list of indices of neighbors for particle i
    double h,                                             // Smoothing length
    double xsph_eps                                       // Epsilon for XSPH
) {
    // Loop over each destination particle 'a' (PySPH d_idx)

        // This corresponds to PySPH's XSPHCorrection.initialize()
        // This accumulates the sum term: -eps * sum_b ( m_b * v_ab / rho_bar_ab * W_ab )
        mysph::vec3<double> xsph_sum_term = {0.0, 0.0, 0.0};

        // Loop over each source particle 'b' (PySPH s_idx) which is a neighbor of 'a'
        for (const auto& pb : neighbors) {
            // Calculate v_ab = v_a - v_b (corresponds to PySPH VIJ)
            mysph::vec3<double> v_ab = pa.v - pb.v;

            // Calculate r_ab = r_a - r_b for the kernel
            mysph::vec3<double> r_ab = pa.r - pb.r;

            // Kernel evaluation W_ab (corresponds to PySPH WIJ)
            double W_ab = mysph::kernel(r_ab, h);

            // Average density term 1/rho_bar_ab (corresponds to PySPH RHOIJ1)
            // PySPH usually uses (rho_i + rho_j)/2 for rho_bar_ij
            if (pa.rho + pb.rho < 1e-9) continue; // Avoid division by zero if densities are very low
            double inv_rho_bar_ab = 2.0 / (pa.rho + pb.rho);

            // The 'tmp' term from PySPH scaled by v_ab
            // tmp = -eps * m_b * W_ab * inv_rho_bar_ab
            // term_contribution = tmp * v_ab
            double factor = -xsph_eps * pb.m * W_ab * inv_rho_bar_ab;
            
            xsph_sum_term[0] += factor * v_ab[0];
            xsph_sum_term[1] += factor * v_ab[1];
            xsph_sum_term[2] += factor * v_ab[2];
        }

        return xsph_sum_term;
}

void compute_artificial_viscosity(
    Particle<double>& pa,                    // All particles
    const std::vector<Particle<double>>& neighbors, // neighbors[i] is list of indices
    double h,
    double visc_alpha,                                   // Viscosity alpha parameter
    double visc_beta,                                    // Viscosity beta parameter
    double eta_factor = 0.1                              // For eta^2 = (eta_factor * h_ab)^2
) {
    const double eta_factor_sq = eta_factor * eta_factor;

    pa.Fv = {0.0, 0.0, 0.0};
    // Loop over each source particle 'b' (PySPH s_idx) which is a neighbor of 'a'
    for (const auto& pb : neighbors) {
        // Relative position: r_ab = r_a - r_b (PySPH XIJ)
        mysph::vec3<double> r_ab = pa.r - pb.r;
        double r_ab_mag_sq = r_ab * r_ab; // PySPH R2IJ
        if (r_ab_mag_sq < 1e-12) continue; // Avoid division by zero if particles are coincident
        // Relative velocity: v_ab = v_a - v_b (PySPH VIJ)
        mysph::vec3<double> v_ab = pa.v - pb.v;
        // v_ab . r_ab (PySPH vijdotxij)
        double v_ab_dot_r_ab = v_ab * r_ab;
        double Pi_ab = 0.0; // Artificial viscosity term
        // Only apply viscosity if particles are approaching (v_ab . r_ab < 0)
        if (v_ab_dot_r_ab < 0.0) {
            // eta_sq parameter (PySPH EPS, often (0.1*h_ab)^2)
            double eta_sq = eta_factor_sq * h * h;
            // phi_ab term (PySPH muij)
            // phi_ab = (h_ab * v_ab . r_ab) / (|r_ab|^2 + eta^2)
            double phi_ab = (h * v_ab_dot_r_ab) / (r_ab_mag_sq + eta_sq);
            // Average speed of sound c_bar_ab (PySPH cij)
            double c_bar_ab = 0.5 * (pa.cs + pb.cs);
            // Average density rho_bar_ab = (rho_a + rho_b) / 2.0
            // PySPH RHOIJ1 = 1.0 / rho_bar_ab
            double rho_bar_ab = 0.5 * (pa.rho + pb.rho);
            if (std::abs(rho_bar_ab) < 1e-9) continue; // Avoid division by zero
            // Calculate Pi_ab
            Pi_ab = pb.m * (-visc_alpha * c_bar_ab * phi_ab + visc_beta * phi_ab * phi_ab) / rho_bar_ab;
        }

        pa.Fv = pa.Fv - Pi_ab * pa.m * mysph::grad_kernel(r_ab, h);
    }
}

