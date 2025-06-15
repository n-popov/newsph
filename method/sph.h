#pragma once

#include "../utils/common.h"
#include "particle.h"

void compute_neighbors(
    Particle<double>& pa,
    const std::vector<Particle<double>*>& particles,
    const config::SPHParameters& sph_params
) {
    for (auto j = 0; j < particles.size(); j++) {
        if (mysph::kernel(pa.r - particles[j]->r, sph_params.h, sph_params.kernel) > 0.) {
            pa.neighbors.push_back(particles[j]);
        }
    }
}

mysph::vec3<double> compute_xsph_corrected_velocities(
    const Particle<double>& pa,                     
    const std::vector<Particle<double>*>& neighbors,
    const config::SPHParameters& sph_params                                     
) {
        mysph::vec3<double> xsph_sum_term = {0.0, 0.0, 0.0};

        for (const auto ppb : neighbors) {
            const auto& pb = *ppb;
            auto v_ab = pa.vstar - pb.vstar;
            auto r_ab = pa.r - pb.r;

            auto W_ab = mysph::kernel(r_ab, sph_params.h, sph_params.kernel);

            if (pa.rho + pb.rho < 1e-9) continue;
            auto inv_rho_bar_ab = 2.0 / (pa.rho + pb.rho);
            auto factor = -sph_params.xsph_eps * pb.m * W_ab * inv_rho_bar_ab;
    
            xsph_sum_term = xsph_sum_term + factor * v_ab;
        }

        return xsph_sum_term;
}

void compute_density(
    Particle<double>& pa,
    const std::vector<Particle<double>*>& neighbors,
    const config::SPHParameters& sph_params
) {
    pa.rho_sph = 0.0;
            
    for (auto ppb : neighbors) {
        pa.rho_sph += ppb->m * mysph::kernel(pa.r - ppb->r, sph_params.h, sph_params.kernel);
    }
}

void compute_correction_factor(
    Particle<double>& pa,
    const std::vector<Particle<double>*>& neighbors,
    const config::SPHParameters& sph_params
) {
    pa.cf = 0.;
            
    for (auto ppb : neighbors) {
        pa.cf += (ppb->m * mysph::kernel(pa.r - ppb->r, sph_params.h, sph_params.kernel)) / ppb->rho;
    }
}

void compute_artificial_viscosity(
    Particle<double>& pa,
    const std::vector<Particle<double>*>& neighbors,
    const config::SPHParameters& sph_params                                  
) {
    auto& eta_factor = sph_params.eta_factor;
    auto& h = sph_params.h;
    const double eta_factor_sq = eta_factor * eta_factor;

    pa.Fv = {0.0, 0.0, 0.0};
    for (const auto ppb : neighbors) {
        const auto& pb = *ppb;
        mysph::vec3<double> r_ab = pa.r - pb.r;
        auto r_ab_mag_sq = r_ab * r_ab;
        if (r_ab_mag_sq < 1e-12) continue;
        mysph::vec3<double> v_ab = pa.v - pb.v;
        auto v_ab_dot_r_ab = v_ab * r_ab;
        auto Pi_ab = 0.0;

        if (v_ab_dot_r_ab < 0.0) {
            auto eta_sq = eta_factor_sq * h * h;
            auto phi_ab = (h * v_ab_dot_r_ab) / (r_ab_mag_sq + eta_sq);
            auto c_bar_ab = 0.5 * (pa.cs + pb.cs);
            auto rho_bar_ab = 0.5 * (pa.rho + pb.rho);
            if (std::abs(rho_bar_ab) < 1e-9) continue;
            Pi_ab = (-sph_params.avisc_alpha * c_bar_ab * phi_ab + sph_params.avisc_beta * phi_ab * phi_ab) / rho_bar_ab;
        }

        pa.a = pa.a - Pi_ab * pb.m * mysph::grad_kernel(r_ab, h, sph_params.kernel);
    }
}



int get_nearest_particle(
    const mysph::Particle<double>& p, const std::vector<mysph::Particle<double>>& neighbors
) {
    auto min_distance = 0.;
    auto nearest = -1;
    for (auto i = 0; i < neighbors.size(); i++) {
        if (neighbors[i].is_fake) continue;

        auto distance = mysph::abs(p.r - neighbors[i].r);
        if (distance < min_distance) {
            nearest = i;
            min_distance = distance;
        }
    }

    return nearest;
}
