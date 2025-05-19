#pragma once

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
    pi.v00 = pi.v01 = pi.v02 = pi.v10 = pi.v11 = pi.v12 = pi.v20 = pi.v21 = pi.v22 = 0.0;

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