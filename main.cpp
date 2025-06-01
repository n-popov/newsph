#include <iostream>
#include <numeric>
#include <vector>
#include <filesystem>
#include <fstream>

#include "utils/vtk.h"
#include "method/kernel.h"
#include "method/particle.h"
#include "utils/common.h"
#include "utils/helpers.h"
#include "method/physics.h"
#include "method/sph.h"

#include "config/simulation_config.h"
#include "geometry/projectile.h"
#include "geometry/plate.h"

int main(int argc, char* argv[]) {
    // Load configuration
    std::string config_file = "config.json";
    if (argc > 1) {
        config_file = argv[1];
    }
    
    config::SimulationConfig config;
    
    try {
        config.load_from_file(config_file);
        std::cout << "Loaded configuration from " << config_file << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error loading configuration: " << e.what() << std::endl;
        return 1;
    }
    
    // Get simulation parameters
    const auto& sim_params = config.simulation_params;
    const auto& sph_params = config.sph_params;
    const auto& gmsh_params = config.gmsh_params;
    
    // Create particles
    std::vector<mysph::Particle<double>> plate_particles;
    std::vector<mysph::Particle<double>> projectile_particles;

    if (gmsh_params.is_enabled) {
        plate_particles = geometry::PlateGenerator::create_plate_from_gmsh(config);
        projectile_particles = geometry::ProjectileGenerator::create_projectile_from_gmsh(config);
    } else {
        plate_particles = geometry::PlateGenerator::create_plate(config);
        projectile_particles = geometry::ProjectileGenerator::create_projectile(config);
    } 

    std::cout << "Created " << plate_particles.size() << " plate particles and " 
              << projectile_particles.size() << " projectile particles\n";
    
    // Combine particles for simulation
    std::vector<mysph::Particle<double>> particles;
    particles.reserve(plate_particles.size() + projectile_particles.size());
    particles.insert(particles.end(), projectile_particles.begin(), projectile_particles.end());
    particles.insert(particles.end(), plate_particles.begin(), plate_particles.end());
    
    int step = 0;
    double time = 0.0;
    
    std::filesystem::create_directory("output");
    
    while (time < sim_params.max_time) {
        std::cout << "Step " << step << ", Time: " << time << " s\n";
        
        auto next_particles = particles;
        
        for (size_t i = 0; i < particles.size(); i++) {
            particles[i].rho = 0.0;
            
            for (size_t j = 0; j < particles.size(); j++) {
                particles[i].rho += particles[j].m * mysph::kernel(particles[i].r - particles[j].r, sph_params.h);
            }
        }

        if (step == 0) {
            std::string vtk_filename = "output/impact-0.vtp";
            write_particles_vtk(vtk_filename, particles);
        }
        
        std::vector<std::vector<mysph::Particle<double>>> neighbors(particles.size());
        for (size_t i = 0; i < particles.size(); i++) {
            for (size_t j = 0; j < particles.size(); j++) {
                if (i != j && mysph::abs(particles[i].r - particles[j].r) <= 2 * sph_params.h) {
                    neighbors[i].push_back(particles[j]);
                }
            }
        }

        for (auto& p : particles) {
            compute_eos_stiffened_gas(p, config);
        }
        
        for (size_t i = 0; i < particles.size(); i++) {
            compute_velocity_gradient(particles[i], neighbors[i], sph_params.h);
        }
        
        for (auto& p : particles) {
            compute_stress_rate_and_artificial_terms(p, sim_params.dt, 0.1);
        }

        for (size_t i = 0; i < particles.size(); i++) {
            compute_artificial_viscosity(particles[i], neighbors[i], sph_params.h, 
                sph_params.avisc_alpha, sph_params.avisc_beta, sph_params.avisc_eta);
        }
        
        for (size_t i = 0; i < particles.size(); i++) {
            compute_force(particles[i], neighbors[i], sph_params);
        }
        
        // forces
        for (size_t i = 0; i < particles.size(); i++) {
            particles[i].v = particles[i].v + (particles[i].F + particles[i].Fv) * (sim_params.dt / particles[i].m);
        }

        // correction
        for (size_t i = 0; i < particles.size(); i++) {
            next_particles[i].v = particles[i].v + compute_xsph_corrected_velocities(particles[i], neighbors[i], sph_params.h, sph_params.xsph_eps);
        }

        // compute next
        for (size_t i = 0; i < particles.size(); i++) {
            next_particles[i].r = particles[i].r + next_particles[i].v * sim_params.dt;
        }

        if (step % sim_params.output_frequency == 0) {
            std::string vtk_filename = "output/impact-" + std::to_string(step + 1) + ".vtp";
            write_particles_vtk(vtk_filename, particles);

            std::string eval_filename = "output/eval-" + std::to_string(step) + ".txt";
            write_full_particle_data(eval_filename, particles, time, step);
        }

        particles = next_particles;
        
        time += sim_params.dt;
        step++;
    }
    
    return 0;
}
