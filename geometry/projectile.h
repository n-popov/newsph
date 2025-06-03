#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <gmsh.h>

#include "../method/particle.h"
#include "../config/simulation_config.h"
#include "../utils/helpers.h"

namespace geometry {

class ProjectileGenerator {
public:
    static std::vector<mysph::Particle<double>> create_projectile(
        const config::SimulationConfig& config) {
        
        std::vector<mysph::Particle<double>> particles;
        
        const auto& proj_params = config.projectile_params;
        const auto& props = config.get_material_properties(proj_params.material);
        double dx = config.sph_params.dx;
        double r_geom = proj_params.radius;
        double v_s = proj_params.velocity;
    

        for (double x = -r_geom; x <= r_geom; x += dx) {
            for (double y = -r_geom; y <= r_geom; y += dx) {
                for (double z = -r_geom; z <= r_geom; z += dx) {
                    double d_sq = x*x + y*y + z*z;
                    if (d_sq <= r_geom*r_geom) {
                        mysph::Particle<double> p;
                        p.r = {x - (r_geom + dx), y + r_geom, z};
                        p.v = {v_s, 0.0, 0.0};
                        p.m = dx * dx * dx * props.density;
                        p.rho0 = props.density;
                        p.rho = props.density;
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
        
        std::cout << "Created " << particles.size() << " projectile particles using grid method\n";
        return particles;
    }

    static std::vector<mysph::Particle<double>> create_projectile_from_gmsh(
        const config::SimulationConfig& config) {
        
        std::vector<mysph::Particle<double>> particles;
        const auto& proj_params = config.projectile_params;
        const auto& props = config.get_material_properties(proj_params.material);
        double dx_mesh = config.sph_params.dx * 2;
        double r_geom = proj_params.radius;
        double v_s = proj_params.velocity;

        double offset_x = -(r_geom + dx_mesh);
        double offset_y = r_geom;
        double offset_z = 0.0;

        gmsh::initialize();
        gmsh::option::setNumber("General.Terminal", config.gmsh_params.verbose ? 1 : 0); // 1 for verbose, 0 for quiet

        try {
            gmsh::model::add("projectile");

            int sphere_tag = gmsh::model::occ::addSphere(0, 0, 0, r_geom);
            gmsh::model::occ::synchronize();

            gmsh::model::addPhysicalGroup(3, {1}, -1, "projectile_volume");

            gmsh::option::setNumber("Mesh.MeshSizeMin", dx_mesh);
            gmsh::option::setNumber("Mesh.MeshSizeMax", dx_mesh);
            gmsh::option::setNumber("Mesh.Algorithm", 6); // Frontal-Delaunay
            gmsh::option::setNumber("Mesh.Algorithm3D", 1); // Delaunay
            gmsh::option::setNumber("Mesh.ElementOrder", 1); // Linear

            gmsh::model::mesh::generate(3); 

            std::vector<double> barycenters;
            gmsh::model::mesh::getBarycenters(4, -1, false, true, barycenters);

            auto size = barycenters.size() / 3;
            auto mass = props.density * 4 * M_PI * r_geom * r_geom * r_geom / (3 * size);

            particles.reserve(size);
            for (size_t i = 0; i < size; ++i) {
                mysph::Particle<double> p;
                p.r = {
                    barycenters[i*3] + offset_x,
                    barycenters[i*3+1] + offset_y,
                    barycenters[i*3+2] + offset_z
                };
                p.v = {v_s, 0.0, 0.0};
                p.m = mass;
                p.rho0 = props.density;
                p.rho = props.density;
                p.cs = props.sound_speed;
                p.material = 1; 
                p.G = props.shear_modulus;
                p.Yo = props.yield_strength;
                p.k = props.bulk_modulus;
                p.nu = props.poisson_ratio;
                particles.push_back(p);
            }

            if (config.gmsh_params.write_mesh_file) {
                gmsh::write("projectile_mesh.msh");
                std::cout << "Wrote projectile_mesh.msh" << std::endl;
            }

        } catch (const std::exception& e) {
            std::cerr << "Gmsh error in projectile generation: " << e.what() << std::endl;
            gmsh::finalize();
            throw e;
        }

        gmsh::finalize();
        
        std::cout << "Created " << particles.size() << " projectile particles using Gmsh\n";
        return particles;
    }
};

} // namespace geometry
