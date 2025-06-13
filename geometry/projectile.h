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
                const auto* nodes = &tet_node_tags[i * 4];
                const auto& v0 = node_coord_map[nodes[0]];
                const auto& v1 = node_coord_map[nodes[1]];
                const auto& v2 = node_coord_map[nodes[2]];
                const auto& v3 = node_coord_map[nodes[3]];

                volumes[i] = std::abs((v1 - v0) * vecmul(v2 - v0, v3 - v0)) / 6.;
            }

            particles.reserve(particles.size() + tet_element_tags.size());
            for (size_t i = 0; i < tet_element_tags.size(); ++i) {
                mysph::Particle<double> p;
                p.r = {
                    barycenters[i*3] + offset_x,
                    barycenters[i*3+1] + offset_y,
                    barycenters[i*3+2] + offset_z
                };
                p.v = {v_s, 0.0, 0.0};
                p.m = volumes[i] * props.density;
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
