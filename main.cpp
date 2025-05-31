#include <iostream>
#include <numeric>
#include <vector>
#include <filesystem>

#include "utils/vtk.h"
#include "utils/kernel.h"
#include "utils/particle.h"
#include "utils/common.h"
#include "utils/helpers.h"
#include "utils/physics.h"
#include "utils/gemini_helpers.h"

double max_time = 2e-6;
double dt = 1e-8;

// Aluminum
const double ro1 = 2700.0;     
const double C1 = 5328.0;      
const double S1 = 1.338;       
const double gamma1 = 2.0;     
const double G1 = 2.76e7;      
const double Yo1 = 0.3e6;      
const double E1 = ro1 * C1 * C1;
const double cs1 = std::sqrt(E1/ro1);

// Steel
const double ro2 = 7900;       
const double C2 = 4600.0;      
const double S2 = 1.490;       
const double gamma2 = 2.0;     
const double G2 = 0.530e7;     
const double Yo2 = 0.979e6;    
const double E2 = ro2 * C2 * C2;
const double cs2 = std::sqrt(E2/ro2);

// SPH parameters
const double hdx = 1.3;
const double dx = 0.0001 * 1.5;     
const double h = dx * hdx;
const double avisc_alpha = 1.0;
const double avisc_beta = 1.5; 
const double avisc_eta = 0.1;  
const double xsph_eps = 0.5;   

// Projectile parameters
const double r = 0.001;        
const double v_s = 2200.0;     

// Function to create projectile particles (sphere) - reduced resolution
std::vector<mysph::Particle<double>> create_projectile() {
    std::vector<mysph::Particle<double>> particles;
    
    // Use a larger particle spacing for testing
    double test_dx = dx;
    
    for (double x = -r; x <= r; x += test_dx) {
        for (double y = -r; y <= r; y += test_dx) {
            for (double z = -r; z <= r; z += test_dx) {
                double d = x*x + y*y + z*z;
                if (d <= r*r) {
                    mysph::Particle<double> p;
                    p.r = {x - (r + test_dx), y + r, z};
                    p.v = {v_s, 0.0, 0.0};
                    p.m = test_dx * test_dx * test_dx * ro2;
                    p.rho = p.rho0 = ro2;
                    p.cs = cs2;
                    p.material = 1; 
                    p.G = G2;
                    p.Yo = Yo2;
                    p.k = E2;
                    p.nu = 0.3; // Poisson ratio 
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
    
    const double plate_length = 0.005;    // 0.025
    const double plate_width = 0.005;     // 0.03
    const double plate_thickness = 0.0005; // 0.0025
    
    double test_dx = dx;
    
    for (double x = 0; x <= plate_thickness; x += test_dx) {
        for (double y = 0; y <= plate_length; y += test_dx) {
            for (double z = -plate_width/2; z <= plate_width/2; z += test_dx) {
                mysph::Particle<double> p;
                p.r = {x, y, z};
                p.v = {0.0, 0.0, 0.0};
                p.m = test_dx * test_dx * test_dx * ro1;
                p.rho = p.rho0 = ro1;
                p.cs = cs1;
                p.material = 0; 
                p.G = G1;
                p.Yo = Yo1;
                p.k = E1;
                p.nu = 0.33; // Poisson ratio 
                particles.push_back(p);
            }
        }
    }
    
    std::cout << "Created " << particles.size() << " plate particles\n";
    return particles;
}

// Stiffened Gas EOS
void compute_pressure(mysph::Particle<double>& p) {
    double gamma, c0, rho0;
    
    if (p.material == 0) { // plate (aluminum)
        gamma = gamma1;
        c0 = C1;
        rho0 = ro1;
    } else { // projectile (steel)
        gamma = gamma2;
        c0 = C2;
        rho0 = ro2;
    }
    
    // Stiffened gas EOS: p = (gamma-1)*rho*e - gamma*p0
    double eta = p.rho / rho0;
    p.p = rho0 * c0 * c0 * (eta - 1) * (1 + (gamma - 1)/2 * (eta - 1));
}

// Hooke's law for deviatoric stress rate
void compute_stress_rate(mysph::Particle<double>& p, double dt) {
    // Calculate deviatoric strain rate
    double div_v = p.v00 + p.v11 + p.v22;
    double eps00 = p.v00 - div_v/3.0;
    double eps01 = 0.5 * (p.v01 + p.v10);
    double eps02 = 0.5 * (p.v02 + p.v20);
    double eps11 = p.v11 - div_v/3.0;
    double eps12 = 0.5 * (p.v12 + p.v21);
    double eps22 = p.v22 - div_v/3.0;
    
    // Update deviatoric stress using Hooke's law
    p.s00 += 2 * p.G * eps00 * dt;
    p.s01 += 2 * p.G * eps01 * dt;
    p.s02 += 2 * p.G * eps02 * dt;
    p.s11 += 2 * p.G * eps11 * dt;
    p.s12 += 2 * p.G * eps12 * dt;
    p.s22 += 2 * p.G * eps22 * dt;
    
    // Apply plasticity
    apply_von_mises_plasticity(p);
}



int main() {
    // Create particles
    auto plate_particles = create_plate();
    // auto plate_particles = std::vector<mysph::Particle<double>>{};
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
            compute_eos_stiffened_gas(p, gamma1, C1, ro1, gamma2, C2, ro2);
        }
        
        for (auto& p : particles) {
            // compute_stress_rate(p, dt);
            compute_stress_rate_and_artificial_terms(p, dt, 0.1);
        }
        
        for (size_t i = 0; i < particles.size(); i++) {
            auto& pi = particles[i];
            
            pi.F = {0.0, 0.0, 0.0};
            pi.Fv = {0.0, 0.0, 0.0};
            
            for (auto& pj : neighbors[i]) {
                auto rij = pi.r - pj.r;
                auto vij = pi.v - pj.v;
                double r = mysph::abs(rij);
                
                // if (r > 0) {
                //     auto eij = rij * (1.0/r);
                //     auto grad_k = mysph::grad_kernel(rij, h);
                    
                //     double dot_prod = vij[0]*rij[0] + vij[1]*rij[1] + vij[2]*rij[2];
                //     double mu_ij = 0.0;
                    
                //     if (dot_prod < 0) {
                //         double c_ij = 0.5 * (pi.cs + pj.cs);
                //         mu_ij = h * dot_prod / (r*r + avisc_eta*h*h);
                //         double pi_ij = (-avisc_alpha * mu_ij * c_ij + avisc_beta * mu_ij * mu_ij) / 
                //                     (0.5 * (pi.rho + pj.rho));
                        
                //         pi.Fv = pi.Fv - grad_k * (pj.m * pi_ij);
                //     }
                    
                //     double stress_term_i[3] = {
                //         (pi.s00 - pi.p) * grad_k[0] + pi.s01 * grad_k[1] + pi.s02 * grad_k[2],
                //         pi.s01 * grad_k[0] + (pi.s11 - pi.p) * grad_k[1] + pi.s12 * grad_k[2],
                //         pi.s02 * grad_k[0] + pi.s12 * grad_k[1] + (pi.s22 - pi.p) * grad_k[2]
                //     };
                    
                //     double stress_term_j[3] = {
                //         (pj.s00 - pj.p) * grad_k[0] + pj.s01 * grad_k[1] + pj.s02 * grad_k[2],
                //         pj.s01 * grad_k[0] + (pj.s11 - pj.p) * grad_k[1] + pj.s12 * grad_k[2],
                //         pj.s02 * grad_k[0] + pj.s12 * grad_k[1] + (pj.s22 - pj.p) * grad_k[2]
                //     };
                    
                //     double density_factor_i = pj.m / (pi.rho * pi.rho);
                //     double density_factor_j = pj.m / (pj.rho * pj.rho);
                    
                //     pi.F[0] -= density_factor_i * stress_term_i[0] + density_factor_j * stress_term_j[0];
                //     pi.F[1] -= density_factor_i * stress_term_i[1] + density_factor_j * stress_term_j[1];
                //     pi.F[2] -= density_factor_i * stress_term_i[2] + density_factor_j * stress_term_j[2];
                // }

                // new impl
                
                double r_dist = mysph::abs(rij);

                if (r_dist > 1e-9 * h) {
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

                    // Now, the force contribution from this neighbor is ma * acc_contrib_vec
                    pi.F[0] += ma * acc_contrib_vec[0];
                    pi.F[1] += ma * acc_contrib_vec[1];
                    pi.F[2] += ma * acc_contrib_vec[2];

                } // end if (r_dist > 0)

                
                
                // auto den = pj.rho * mysph::abs(rij);
                // if (den == 0) continue;
                // pi.Fv = pi.Fv - (pi.v - pj.v) * (pj.m * 2 * mysph::abs(mysph::grad_kernel(rij, h)) / (den));
                // pi.F = pi.F - mysph::grad_kernel(rij, h) * (pj.m * (pi.p / std::pow(pi.rho, 2) + pj.p / std::pow(pj.rho, 2)));
            }
        }
        
        for (size_t i = 0; i < particles.size(); i++) {
            // next_particles[i].v = particles[i].v + (particles[i].Fv) * (dt / particles[i].m); 
            // next_particles[i].v = next_particles[i].v + particles[i].F * (dt / particles[i].m);
            next_particles[i].v = particles[i].v + particles[i].F * (dt / particles[i].m);
            
            next_particles[i].r = particles[i].r + next_particles[i].v * dt;
        }
        
        // XSPH 
        for (size_t i = 0; i < particles.size(); i++) {
            mysph::vec3<double> xsph_correction = {0.0, 0.0, 0.0};
            
            for (auto& pj : neighbors[i]) {
                auto rij = particles[i].r - pj.r;
                auto vij = pj.v - particles[i].v;
                
                double w = mysph::kernel(rij, h);
                double weight = 2.0 * pj.m / (particles[i].rho + pj.rho) * w;
                
                xsph_correction[0] += weight * vij[0];
                xsph_correction[1] += weight * vij[1];
                xsph_correction[2] += weight * vij[2];
            }
            
            next_particles[i].r[0] += xsph_eps * xsph_correction[0] * dt;
            next_particles[i].r[1] += xsph_eps * xsph_correction[1] * dt;
            next_particles[i].r[2] += xsph_eps * xsph_correction[2] * dt;
        }
        
        particles = next_particles;

        if (step % 1 == 0) {
            std::string vtk_filename = "output/impact-" + std::to_string(step + 1) + ".vtp";
            write_particles_vtk(vtk_filename, particles);

            std::string eval_filename = "output/eval-" + std::to_string(step) + ".txt";
            write_full_particle_data(eval_filename, particles, time, step);
        }

        
        time += dt;
        step++;
    }
    
    return 0;
}
