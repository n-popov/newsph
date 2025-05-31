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
