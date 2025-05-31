#include <iostream>
#include <numeric>
#include <vector>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

#include "utils/vtk.h"
#include "utils/kernel.h"
#include "utils/particle.h"
#include "utils/common.h"
#include "utils/helpers.h"
#include "utils/physics.h"
#include "utils/sph.h"

using json = nlohmann::json;

// Global configuration object
json config;

// Material properties structure
struct MaterialProperties {
    double density;
    double sound_speed_coefficient;
    double hugoniot_slope;
    double gruneisen_gamma;
    double shear_modulus;
    double yield_strength;
    double poisson_ratio;
    double bulk_modulus;
    double sound_speed;
};

// Global material properties
MaterialProperties aluminum_props;
MaterialProperties steel_props;

// Load configuration from JSON file
void load_config(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open config file: " + filename);
    }
    file >> config;
    file.close();
    
    // Parse material properties
    auto& al = config["materials"]["aluminum"];
    aluminum_props.density = al["density"];
    aluminum_props.sound_speed_coefficient = al["sound_speed_coefficient"];
    aluminum_props.hugoniot_slope = al["hugoniot_slope"];
    aluminum_props.gruneisen_gamma = al["gruneisen_gamma"];
    aluminum_props.shear_modulus = al["shear_modulus"];
    aluminum_props.yield_strength = al["yield_strength"];
    aluminum_props.poisson_ratio = al["poisson_ratio"];
    aluminum_props.bulk_modulus = aluminum_props.density * 
        aluminum_props.sound_speed_coefficient * aluminum_props.sound_speed_coefficient;
    aluminum_props.sound_speed = std::sqrt(aluminum_props.bulk_modulus / aluminum_props.density);
    
    auto& st = config["materials"]["steel"];
    steel_props.density = st["density"];
    steel_props.sound_speed_coefficient = st["sound_speed_coefficient"];
    steel_props.hugoniot_slope = st["hugoniot_slope"];
    steel_props.gruneisen_gamma = st["gruneisen_gamma"];
    steel_props.shear_modulus = st["shear_modulus"];
    steel_props.yield_strength = st["yield_strength"];
    steel_props.poisson_ratio = st["poisson_ratio"];
    steel_props.bulk_modulus = steel_props.density * 
        steel_props.sound_speed_coefficient * steel_props.sound_speed_coefficient;
    steel_props.sound_speed = std::sqrt(steel_props.bulk_modulus / steel_props.density);
}

// Function to create projectile particles (sphere)
std::vector<mysph::Particle<double>> create_projectile() {
    std::vector<mysph::Particle<double>> particles;
    
    double r = config["projectile"]["radius"];
    double v_s = config["projectile"]["velocity"];
    double dx = config["sph_parameters"]["dx"];
    
    const MaterialProperties& props = (config["projectile"]["material"] == "steel") ? 
        steel_props : aluminum_props;
    
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

std::vector<mysph::Particle<double>> create_plate() {
    std::vector<mysph::Particle<double>> particles;
    
    double plate_length = config["plate"]["length"];
    double plate_width = config["plate"]["width"];
    double plate_thickness = config["plate"]["thickness"];
    double dx = config["sph_parameters"]["dx"];
    
    const MaterialProperties& props = (config["plate"]["material"] == "aluminum") ? 
        aluminum_props : steel_props;
    
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



int main(int argc, char* argv[]) {
    // Load configuration
    std::string config_file = "config.json";
    if (argc > 1) {
        config_file = argv[1];
    }
    
    try {
        load_config(config_file);
        std::cout << "Loaded configuration from " << config_file << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error loading configuration: " << e.what() << std::endl;
        return 1;
    }
    
    // Get simulation parameters
    double max_time = config["simulation"]["max_time"];
    double dt = config["simulation"]["dt"];
    int output_frequency = config["simulation"]["output_frequency"];
    
    // Get SPH parameters
    double hdx = config["sph_parameters"]["hdx"];
    double dx = config["sph_parameters"]["dx"];
    double h = dx * hdx;
    double avisc_alpha = config["sph_parameters"]["avisc_alpha"];
    double avisc_beta = config["sph_parameters"]["avisc_beta"];
    double avisc_eta = config["sph_parameters"]["avisc_eta"];
    double xsph_eps = config["sph_parameters"]["xsph_eps"];
    
    // Create particles
    auto plate_particles = create_plate();
    auto projectile_particles = create_projectile();
    
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
    
    while (time < max_time) {
        std::cout << "Step " << step << ", Time: " << time << " s\n";
        
        auto next_particles = particles;
        
        for (size_t i = 0; i < particles.size(); i++) {
            particles[i].rho = 0.0;
            
            for (size_t j = 0; j < particles.size(); j++) {
                particles[i].rho += particles[j].m * mysph::kernel(particles[i].r - particles[j].r, h);
            }
        }

        if (step == 0) {
            std::string vtk_filename = "output/impact-0.vtp";
            write_particles_vtk(vtk_filename, particles);
        }
        
        std::vector<std::vector<mysph::Particle<double>>> neighbors(particles.size());
        for (size_t i = 0; i < particles.size(); i++) {
            for (size_t j = 0; j < particles.size(); j++) {
                if (i != j && mysph::abs(particles[i].r - particles[j].r) <= 2 * h) {
                    neighbors[i].push_back(particles[j]);
                }
            }
        }
        
        for (size_t i = 0; i < particles.size(); i++) {
            compute_velocity_gradient(particles, i, neighbors[i], h);
        }
        
        for (auto& p : particles) {
            compute_eos_stiffened_gas(p, aluminum_props.gruneisen_gamma, 
                aluminum_props.sound_speed_coefficient, aluminum_props.density,
                steel_props.gruneisen_gamma, steel_props.sound_speed_coefficient, 
                steel_props.density);
        }
        
        for (auto& p : particles) {
            compute_stress_rate_and_artificial_terms(p, dt, 0.1);
        }

        for (size_t i = 0; i < particles.size(); i++) {
            compute_artificial_viscosity(particles[i], neighbors[i], h, avisc_alpha, avisc_beta, avisc_eta);
        }
        
        for (size_t i = 0; i < particles.size(); i++) {
            auto& pi = particles[i];
            
            pi.F = {0.0, 0.0, 0.0};
            
            for (auto& pj : neighbors[i]) {
                auto rij = pi.r - pj.r;
                auto vij = pi.v - pj.v;
                double r = mysph::abs(rij);

                if (r > 1e-9 * h) {
                    mysph::vector3d DWIJ = mysph::grad_kernel(rij, h);
                    double WIJ = mysph::kernel(rij, h);

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

                    auto wdeltap = mysph::kernel(mysph::vector3d{r, r, r}, h);
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
            particles[i].v = particles[i].v + (particles[i].F + particles[i].Fv) * (dt / particles[i].m);
        }

        // correction
        for (size_t i = 0; i < particles.size(); i++) {
            next_particles[i].v = particles[i].v + compute_xsph_corrected_velocities(particles[i], neighbors[i], h, xsph_eps);
        }

        // compute next
        for (size_t i = 0; i < particles.size(); i++) {
            next_particles[i].r = particles[i].r + next_particles[i].v * dt;
        }

        if (step % output_frequency == 0) {
            std::string vtk_filename = "output/impact-" + std::to_string(step + 1) + ".vtp";
            write_particles_vtk(vtk_filename, particles);

            std::string eval_filename = "output/eval-" + std::to_string(step) + ".txt";
            write_full_particle_data(eval_filename, particles, time, step);
        }

        particles = next_particles;
        
        time += dt;
        step++;
    }
    
    return 0;
}
