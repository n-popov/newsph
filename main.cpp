#include <iostream>
#include <numeric>
#include <vector>
#include <filesystem>
#include <format>
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
        if (sim_params.use_plate){
            plate_particles = geometry::PlateGenerator::create_plate_from_gmsh(config);
        }
        if (sim_params.use_projectile) {
            projectile_particles = geometry::ProjectileGenerator::create_projectile_from_gmsh(config);
        }
    } else {
        if (sim_params.use_plate) {
            plate_particles = geometry::PlateGenerator::create_plate(config);
        }

        if (sim_params.use_projectile) {
            projectile_particles = geometry::ProjectileGenerator::create_projectile(config);
        }
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

    auto max_cf = 0.;
    auto v_max_abs = 0.;

    for (const auto& p: particles) {
        v_max_abs = std::max(v_max_abs, abs(p.v));
    }
    
    while (time < sim_params.max_time) {
        std::cout << "Step " << step << ", Time: " << time << " s\n";
        
        auto next_particles = particles;

        std::vector<mysph::Particle<double>*> particles_pointers(std::size(particles));
        std::transform(particles.begin(), particles.end(), particles_pointers.begin(), [](auto& particle){return &particle;});

        for (auto& p: particles) {
            p.external_neighbors = &particles_pointers;
            p.neighbors = {};
        }

        parallelize(sim_params.parallelize, compute_neighbors, particles, sph_params);

        parallelize(sim_params.parallelize, compute_correction_factor, particles, sph_params);

        parallelize(sim_params.parallelize, continuity_equation, particles, sph_params);

        if (step == 0) {
            std::string vtk_filename = "output/impact-0.vtp";
            write_particles_vtk(vtk_filename, particles);

            for (auto& p: particles) {
                max_cf = std::max(max_cf, p.cf);
            }
        } else {
            parallelize(sim_params.parallelize, compute_density, particles, sph_params);

            // density selection
            for (auto& p: particles) {
                if (std::abs(p.rho_sph - p.rho) < std::abs(p.rho_cont - p.rho)) {
                    p.rho = p.rho_sph;
                } else {
                    p.rho = p.rho_cont;
                }
            }
        }

        for (auto& p : particles) {
            compute_eos_stiffened_gas(p, config);
        }

        parallelize(sim_params.parallelize, compute_velocity_gradient, particles, sph_params);
        
        for (auto& p : particles) {
            compute_hookes_deviatoric_stress_rate(p);

            compute_monaghan_artificial_stress(p, 0.1);
        }

        parallelize(sim_params.parallelize, compute_artificial_viscosity, particles, sph_params);

        parallelize(sim_params.parallelize, compute_force, particles, sph_params);

        parallelize(sim_params.parallelize, compute_energy_with_stress, particles, sph_params);
        
        // forces
        for (auto i = 0; i < particles.size(); i++) {
            particles[i].vstar = particles[i].v + particles[i].a * sim_params.dt;
        }

        // correction
        for (auto i = 0; i < particles.size(); i++) {
            // high-speed correction
            if (abs(particles[i].vstar) > sph_params.hs_factor * v_max_abs) {
                particles[i].vstar = particles[i].vstar * (v_max_abs / abs(particles[i].vstar));
            }
            particles[i].v = particles[i].vstar + compute_xsph_corrected_velocities(particles[i], particles[i].neighbors, sph_params);
        }

        // compute next
        for (auto i = 0; i < particles.size(); i++) {
            next_particles[i].r = particles[i].r + particles[i].v * sim_params.dt;
            next_particles[i].e = particles[i].e + particles[i].ae * sim_params.dt;
            next_particles[i].stress = particles[i].stress + particles[i].acc_stress * sim_params.dt;
            next_particles[i].rho_cont = particles[i].rho + particles[i].arho * sim_params.dt;

            next_particles[i].rho = particles[i].rho;
        } 

        if (step % sim_params.output_frequency == 0) {
            std::string vtk_filename = "output/impact-" + std::to_string(step + 1) + ".vtp";
            write_particles_vtk(vtk_filename, particles);

            if (sim_params.write_full_evaluation_data) {
                std::string eval_filename = "output/eval-" + std::to_string(step) + ".txt";
                write_full_particle_data(eval_filename, particles, time, step);
            }

            if (sim_params.write_integral_data) {
                std::string eval_filename = "output/integral-" + std::to_string(step) + ".txt";
                write_integral_characteristics(eval_filename, particles, time, step);
            }
        }

        particles = next_particles;
        
        time += sim_params.dt;
        step++;
    }
    
    return 0;
}
