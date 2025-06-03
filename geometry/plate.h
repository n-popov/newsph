#pragma once

#include <vector>
#include <iostream>
#include <stdexcept> 
#include <gmsh.h> 

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
        
        double plate_length = plate_params.length;       // Y   
        double plate_width = plate_params.width;         // Z
        double plate_thickness = plate_params.thickness; // X
        
        for (double x = 0; x <= plate_thickness - dx/2.0; x += dx) { 
            for (double y = 0; y <= plate_length - dx/2.0; y += dx) {
                for (double z = -plate_width/2.0; z <= plate_width/2.0 - dx/2.0; z += dx) {
                    mysph::Particle<double> p;
                    p.r = {x + dx/2.0, y + dx/2.0, z + dx/2.0}; 
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
        
        std::cout << "Created " << particles.size() << " plate particles using grid method\n";
        return particles;
    }

    static std::vector<mysph::Particle<double>> create_plate_from_gmsh(
        const config::SimulationConfig& config) {
        
        std::vector<mysph::Particle<double>> particles;
        const auto& plate_params = config.plate_params;
        const auto& props = config.get_material_properties(plate_params.material);
        double dx_mesh = config.sph_params.dx * 2;
        
        double plate_length = plate_params.length;       // Y
        double plate_width = plate_params.width;         // Z
        double plate_thickness = plate_params.thickness; // X

        double x_start = 0.0;
        double y_start = 0.0;
        double z_start = -plate_width / 2.0;

        gmsh::initialize();
        gmsh::option::setNumber("General.Terminal", config.gmsh_params.verbose ? 1 : 0); // 1 for verbose, 0 for quite

        try {
            gmsh::model::add("plate");

            int box_tag = gmsh::model::occ::addBox(x_start, y_start, z_start, 
                                                   plate_thickness, plate_length, plate_width);
            gmsh::model::occ::synchronize();

            gmsh::model::addPhysicalGroup(3, {1}, -1, "plate_volume");

            gmsh::option::setNumber("Mesh.MeshSizeMin", dx_mesh);
            gmsh::option::setNumber("Mesh.MeshSizeMax", dx_mesh);
            gmsh::option::setNumber("Mesh.Algorithm", 6); 
            gmsh::option::setNumber("Mesh.Algorithm3D", 1); 
            gmsh::option::setNumber("Mesh.ElementOrder", 1);

            gmsh::model::mesh::generate(3);

            std::vector<double> barycenters;
            gmsh::model::mesh::getBarycenters(4, -1, false, true, barycenters);

            auto size = barycenters.size() / 3;
            auto mass = props.density * plate_length * plate_width * plate_thickness / size;

            particles.reserve(size);
            for (size_t i = 0; i < size; ++i) {
                mysph::Particle<double> p;
                    p.r = {barycenters[i*3 + 0], barycenters[i*3 + 1], barycenters[i*3 + 2]};
                    p.v = {0.0, 0.0, 0.0};
                    p.m = mass;
                    p.rho = props.density;
                    p.rho0 = props.density;
                    p.cs = props.sound_speed;
                    p.material = 0; 
                    p.G = props.shear_modulus;
                    p.Yo = props.yield_strength;
                    p.k = props.bulk_modulus;
                    p.nu = props.poisson_ratio;
                    particles.push_back(p);
            }

            if (config.gmsh_params.write_mesh_file) { 
                gmsh::write("plate_mesh.msh");
                std::cout << "Wrote plate_mesh.msh" << std::endl;
            }

        } catch (const std::exception& e) {
            std::cerr << "Gmsh error in plate generation: " << e.what() << std::endl;
            gmsh::finalize();
            throw; 
        }

        gmsh::finalize();
        
        std::cout << "Created " << particles.size() << " plate particles using Gmsh\n";
        return particles;
    }
};

} // namespace geometry
