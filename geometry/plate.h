#ifndef PLATE_GENERATOR_H
#define PLATE_GENERATOR_H

#include <vector>
#include <iostream>
#include "../method/particle.h"
#include "../config/simulation_config.h"

namespace geometry {

class PlateGenerator {
public:
    static std::vector<mysph::Particle<double>> create_plate(
        const config::SimulationConfig& config) {
        
        std::vector<mysph::Particle<double>> particles;
        
        const auto& plate_params = config.plate_params;
        const auto& props = config.get_material_properties(plate_params.material);
        double dx = config.sph_params.dx;
        
        double plate_length = plate_params.length;
        double plate_width = plate_params.width;
        double plate_thickness = plate_params.thickness;
        
        for (double x = 0; x <= plate_thickness; x += dx) {
            for (double y = 0; y <= plate_length; y += dx) {
                for (double z = -plate_width/2; z <= plate_width/2; z += dx) {
                    mysph::Particle<double> p;
                    p.r = {x, y, z};
                    p.v = {0.0, 0.0, 0.0};
                    p.m = dx * dx * dx * props.density;
                    p.rho = p.rho0 = props.density;
                    p.cs = props.sound_speed;
                    p.material = 0; 
                    p.G = props.shear_modulus;
                    p.Yo = props.yield_strength;
                    p.k = props.bulk_modulus;
                    p.nu = props.poisson_ratio;
                    particles.push_back(p);
                }
            }
        }
        
        std::cout << "Created " << particles.size() << " plate particles\n";
        return particles;
    }
};

} // namespace geometry

#endif // PLATE_GENERATOR_H
