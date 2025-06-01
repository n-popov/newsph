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

            std::vector<int> element_types;
            gmsh::model::mesh::getElementTypes(element_types, 3); 

            std::vector<std::size_t> all_node_tags;
            std::vector<double> all_node_coords;
            std::vector<double> parametric_coords;
            gmsh::model::mesh::getNodes(all_node_tags, all_node_coords, parametric_coords, -1, -1, false, false);

            std::unordered_map<std::size_t, std::array<double, 3>> node_coord_map;
            for (size_t i = 0; i < all_node_tags.size(); ++i) {
                node_coord_map[all_node_tags[i]] = {
                    all_node_coords[i*3], all_node_coords[i*3+1], all_node_coords[i*3+2]
                };
            }
            
            std::vector<std::size_t> tet_element_tags;
            std::vector<std::size_t> tet_node_tags;
            gmsh::model::mesh::getElementsByType(4, tet_element_tags, tet_node_tags, -1);

            std::vector<double> volumes(tet_element_tags.size());
            std::vector<double> barycenters;
            gmsh::model::mesh::getBarycenters(4, -1, false, true, barycenters);

            for (size_t i = 0; i < tet_element_tags.size(); ++i) {
                const auto* nodes = &tet_node_tags[i*4];
                const auto& v0 = node_coord_map[nodes[0]];
                const auto& v1 = node_coord_map[nodes[1]];
                const auto& v2 = node_coord_map[nodes[2]];
                const auto& v3 = node_coord_map[nodes[3]];

                const double ax = v1[0]-v0[0], ay = v1[1]-v0[1], az = v1[2]-v0[2];
                const double bx = v2[0]-v0[0], by = v2[1]-v0[1], bz = v2[2]-v0[2];
                const double cx = v3[0]-v0[0], cy = v3[1]-v0[1], cz = v3[2]-v0[2];

                volumes[i] = std::abs(ax*(by*cz-bz*cy) + ay*(bz*cx-bx*cz) + az*(bx*cy-by*cx)) / 6.0;
            }

            particles.reserve(particles.size() + tet_element_tags.size());
            for (size_t i = 0; i < tet_element_tags.size(); ++i) {
                mysph::Particle<double> p;
                    p.r = {barycenters[i*3 + 0], barycenters[i*3 + 1], barycenters[i*3 + 2]};
                    p.v = {0.0, 0.0, 0.0};
                    p.m = volumes[i] * props.density;
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
