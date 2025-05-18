void compute_eos_stiffened_gas(
    mysph::Particle<double>& particle,
    double gamma1, double C1, double ro1, // Material 0 properties
    double gamma2, double C2, double ro2  // Material 1 properties
) {
    double gamma_eos, c0_eos, rho0_eos;

    if (particle.material == 0) { // e.g., plate (aluminum)
        gamma_eos = gamma1;
        c0_eos    = C1;
        rho0_eos  = ro1;
    } else { // e.g., projectile (steel)
        gamma_eos = gamma2;
        c0_eos    = C2;
        rho0_eos  = ro2;
    }

    // ---- Pressure Calculation (Matches PySPH Eq. 17) ----
    // d_p[d_idx] = self.c0*self.c0 * (d_rho[d_idx] - self.r0) + \
    //               (self.gamma - 1.0)*d_rho[d_idx]*d_e[d_idx]
    particle.p = c0_eos * c0_eos * (particle.rho - rho0_eos) + \
                 (gamma_eos - 1.0) * particle.rho * particle.e;

    // ---- Speed of Sound Calculation (Matches PySPH Eq. 21) ----
    // d_cs[d_idx] = sqrt(
    //     self.c0*self.c0 + (self.gamma - 1.0) *
    //     (d_e[d_idx] + d_p[d_idx]/d_rho[d_idx])
    // )
    // Ensure density is positive to avoid division by zero or sqrt of negative.
    if (particle.rho > 1e-9) { // Check for positive density
        double term_in_sqrt = c0_eos * c0_eos +
                              (gamma_eos - 1.0) * (particle.e + particle.p / particle.rho);
        if (term_in_sqrt >= 0) {
            particle.cs = std::sqrt(term_in_sqrt);
        } else {
            // Handle cases where term_in_sqrt might be negative (e.g., due to numerical issues or extreme states)
            // This could indicate a problem or require a specific physical interpretation.
            // For simplicity, setting to c0, but proper handling depends on the simulation context.
            particle.cs = c0_eos;
            // Or, you might want to throw an error or log a warning:
            // std::cerr << "Warning: Negative term in sqrt for cs calculation. Particle density: " << particle.rho << ", energy: " << particle.e << ", pressure: " << particle.p << std::endl;
        }
    } else {
        // Handle zero or negative density if it can occur
        particle.cs = c0_eos; // A fallback, or handle as an error
    }
}

// Assuming 'h' (smoothing length) is available, e.g., passed as argument
void compute_velocity_gradient(std::vector<mysph::Particle<double>>& particles, int i,
                               const std::vector<mysph::Particle<double>>& neighbors,
                               double h) { // Added 'h' as an argument
    auto& pi = particles[i];
    pi.v00 = pi.v01 = pi.v02 = pi.v10 = pi.v11 = pi.v12 = pi.v20 = pi.v21 = pi.v22 = 0.0;

    // Use const auto& for pj if it's not modified
    for (const auto& pj : neighbors) {
        auto rij = pi.r - pj.r;
        auto vij_corrected = pj.v - pi.v; // Changed to pj.v - pi.v to match formula (v_b - v_a)
        auto grad_k = mysph::grad_kernel(rij, h); // grad_k should be ∇_a W_ab

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
