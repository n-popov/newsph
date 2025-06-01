#pragma once

#include <vector>
#include <iostream>
#include "../method/particle.h"
#include "../config/simulation_config.h"

namespace geometry {

class ProjectileGenerator {
public:
    static std::vector<mysph::Particle<double>> create_projectile(
        const config::SimulationConfig& config) {
        
        std::vector<mysph::Particle<double>> particles;
        
        const auto& proj_params = config.projectile_params;
        const auto& props = config.get_material_properties(proj_params.material);
        double dx = config.sph_params.dx;
        double r = proj_params.radius;
        double v_s = proj_params.velocity;
        
        for (double x = -r; x <= r; x += dx) {
            for (double y = -r; y <= r; y += dx) {
                for (double z = -r; z <= r; z += dx) {
                    double d = x*x + y*y + z*z;
                    if (d <= r*r) {
                        mysph::Particle<double> p;
                        p.r = {x - (r + dx), y + r, z};
                        p.v = {v_s, 0.0, 0.0};
                        p.m = dx * dx * dx * props.density;
                        p.rho = p.rho0 = props.density;
                        p.cs = props.sound_speed;
                        p.material = 1; 
                        p.G = props.shear_modulus;
                        p.Yo = props.yield_strength;
                        p.k = props.bulk_modulus;
                        p.nu = props.poisson_ratio;
                        particles.push_back(p);
                    }
                }
            }
        }
        
        std::cout << "Created " << particles.size() << " projectile particles\n";
        return particles;
    }
};

} // namespace geometry
