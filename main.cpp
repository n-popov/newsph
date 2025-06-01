#include <iostream>
#include <numeric>
#include <vector>
#include <filesystem>
#include <fstream>

#include "utils/vtk.h"
#include "utils/kernel.h"
#include "utils/particle.h"
#include "utils/common.h"
#include "utils/helpers.h"
#include "utils/physics.h"
#include "utils/sph.h"

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
    
    // Create particles
    auto plate_particles = geometry::PlateGenerator::create_plate(config);
    auto projectile_particles = geometry::ProjectileGenerator::create_projectile(config);
    
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
        
        for (size_t i = 0; i < particles.size(); i++) {
            compute_velocity_gradient(particles, i, neighbors[i], sph_params.h);
        }
        
        for (auto& p : particles) {
            compute_eos_stiffened_gas(p, 
                config.aluminum_props.gruneisen_gamma, 
                config.aluminum_props.sound_speed_coefficient, 
                config.aluminum_props.density,
                config.steel_props.gruneisen_gamma, 
                config.steel_props.sound_speed_coefficient, 
                config.steel_props.density);
        }
        
        for (auto& p : particles) {
            compute_stress_rate_and_artificial_terms(p, sim_params.dt, 0.1);
        }

        for (size_t i = 0; i < particles.size(); i++) {
            compute_artificial_viscosity(particles[i], neighbors[i], sph_params.h, 
                sph_params.avisc_alpha, sph_params.avisc_beta, sph_params.avisc_eta);
        }
        
        for (size_t i = 0; i < particles.size(); i++) {
            auto& pi = particles[i];
            
            pi.F = {0.0, 0.0, 0.0};
            
            for (auto& pj : neighbors[i]) {
                auto rij = pi.r - pj.r;
                auto vij = pi.v - pj.v;
                double r = mysph::abs(rij);

                if (r > 1e-9 * sph_params.h) {
                    mysph::vector3d DWIJ = mysph::grad_kernel(rij, sph_params.h);
                    double WIJ = mysph::kernel(rij, sph_params.h);

                    double pa = pi.p;
                    double pb = pj.p;
                    double rhoa = pi.rho;
                    double rhob = pj.rho;

                    double rhoa21 = 1.0 / (rhoa * rhoa);
                    double rhob21 = 1.0 / (rhob * rhob);

                    double s00a_dev = pi.s00; double s01a_dev = pi.s01; double s02a_dev = pi.s02;
                    double s11a_dev = pi.s11; double s12a_dev = pi.s12;
                    double s22a_dev = pi.s22;

                    double s00b_dev = pj.s00; double s01b_dev = pj.s01; double s02b_dev = pj.s02;
                    double s11b_dev = pj.s11; double s12b_dev = pj.s12;
                    double s22b_dev = pj.s22;

                    double sigma00a = s00a_dev - pa; double sigma01a = s01a_dev;       double sigma02a = s02a_dev;
                    double sigma10a = s01a_dev;       double sigma11a = s11a_dev - pa; double sigma12a = s12a_dev;
                    double sigma20a = s02a_dev;       double sigma21a = s12a_dev;       double sigma22a = s22a_dev - pa;

                    double sigma00b = s00b_dev - pb; double sigma01b = s01b_dev;       double sigma02b = s02b_dev;
                    double sigma10b = s01b_dev;       double sigma11b = s11b_dev - pb; double sigma12b = s12b_dev;
                    double sigma20b = s02b_dev;       double sigma21b = s12b_dev;       double sigma22b = s22b_dev - pb;

                    double r00_ab = pi.as00 + pj.as00; double r01_ab = pi.as01 + pj.as01; double r02_ab = pi.as02 + pj.as02;
                    double r11_ab = pi.as11 + pj.as11; double r12_ab = pi.as12 + pj.as12;
                    double r22_ab = pi.as22 + pj.as22;

                    double fab_pow_n = 0.0;
                    double art_stress00 = 0.0, art_stress01 = 0.0, art_stress02 = 0.0;
                    double art_stress11 = 0.0, art_stress12 = 0.0;
                    double art_stress22 = 0.0;

                    auto wdeltap = mysph::kernel(mysph::vector3d{r, r, r}, sph_params.h);
                    auto n_art_stress = 2;

                    if (wdeltap > 1e-9) {
                        double fab = WIJ / wdeltap;
                        fab_pow_n = std::pow(fab, n_art_stress);

                        art_stress00 = fab_pow_n * r00_ab;
                        art_stress01 = fab_pow_n * r01_ab;
                        art_stress02 = fab_pow_n * r02_ab;
                        art_stress11 = fab_pow_n * r11_ab;
                        art_stress12 = fab_pow_n * r12_ab;
                        art_stress22 = fab_pow_n * r22_ab;
                    }

                    double mb = pj.m;
                    double ma = pi.m; 
                    mysph::vector3d acc_contrib_vec;

                    double term_xx = (sigma00a * rhoa21 + sigma00b * rhob21 + art_stress00);
                    double term_xy = (sigma01a * rhoa21 + sigma01b * rhob21 + art_stress01);
                    double term_xz = (sigma02a * rhoa21 + sigma02b * rhob21 + art_stress02);
                    acc_contrib_vec[0] = mb * (term_xx * DWIJ[0] + term_xy * DWIJ[1] + term_xz * DWIJ[2]);

                    double term_yx = (sigma10a * rhoa21 + sigma10b * rhob21 + art_stress01);
                    double term_yy = (sigma11a * rhoa21 + sigma11b * rhob21 + art_stress11);
                    double term_yz = (sigma12a * rhoa21 + sigma12b * rhob21 + art_stress12);
                    acc_contrib_vec[1] = mb * (term_yx * DWIJ[0] + term_yy * DWIJ[1] + term_yz * DWIJ[2]);

                    double term_zx = (sigma20a * rhoa21 + sigma20b * rhob21 + art_stress02);
                    double term_zy = (sigma21a * rhoa21 + sigma21b * rhob21 + art_stress12);
                    double term_zz = (sigma22a * rhoa21 + sigma22b * rhob21 + art_stress22);
                    acc_contrib_vec[2] = mb * (term_zx * DWIJ[0] + term_zy * DWIJ[1] + term_zz * DWIJ[2]);

                    pi.F[0] += ma * acc_contrib_vec[0];
                    pi.F[1] += ma * acc_contrib_vec[1];
                    pi.F[2] += ma * acc_contrib_vec[2];
                }
            }
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
