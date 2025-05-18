#include <iostream>
#include <numeric>
#include <vector>
#include <filesystem>

#include "utils/vtk.h"
#include "utils/kernel.h"
#include "utils/particle.h"
#include "utils/common.h"

// Simulation parameters
double max_time = 1e-6;
double dt = 1e-8;

// Material properties - Aluminum (plate)
const double ro1 = 2700.0;     // reference density
const double C1 = 5328.0;      // reference sound speed
const double S1 = 1.338;       // Particle-shock velocity slope
const double gamma1 = 2.0;     // Gruneisen parameter
const double G1 = 2.76e7;      // Shear Modulus
const double Yo1 = 0.3e6;      // Yield stress
const double E1 = ro1 * C1 * C1;
const double cs1 = std::sqrt(E1/ro1);

// Material properties - Steel (projectile)
const double ro2 = 7900;       // reference density
const double C2 = 4600.0;      // reference sound speed
const double S2 = 1.490;       // Particle-shock velocity slope
const double gamma2 = 2.0;     // Gruneisen parameter
const double G2 = 0.530e7;     // Shear Modulus
const double Yo2 = 0.979e6;    // Yield stress
const double E2 = ro2 * C2 * C2;
const double cs2 = std::sqrt(E2/ro2);

// SPH parameters
const double h = 0.00013;      // smoothing length
const double hdx = 1.3;        // h/dx ratio
const double dx = h / hdx;     // particle spacing
const double avisc_alpha = 1.0;// artificial viscosity parameter
const double avisc_beta = 1.5; // artificial viscosity parameter
const double avisc_eta = 0.1;  // artificial viscosity small denominator fix
const double xsph_eps = 0.5;   // XSPH correction strength

// Projectile parameters
const double r = 0.001;        // projectile radius
const double v_s = 2200.0;     // projectile velocity

// Function to create projectile particles (sphere) - reduced resolution
std::vector<mysph::Particle<double>> create_projectile() {
    std::vector<mysph::Particle<double>> particles;
    
    // Use a larger particle spacing for testing
    double test_dx = dx * 2.5;
    
    for (double x = -r; x <= r; x += test_dx) {
        for (double y = -r; y <= r; y += test_dx) {
            for (double z = -r; z <= r; z += test_dx) {
                double d = x*x + y*y + z*z;
                if (d <= r*r) {
                    mysph::Particle<double> p;
                    p.r = {x - (r + test_dx), y - r, z};
                    p.v = {v_s, 2*v_s, 0.0};
                    p.m = test_dx * test_dx * test_dx * ro2;
                    p.rho = p.rho0 = ro2;
                    p.cs = cs2;
                    p.material = 1; // projectile
                    p.G = G2;
                    p.Yo = Yo2;
                    p.k = E2;
                    p.nu = 0.3; // Poisson ratio for steel
                    particles.push_back(p);
                }
            }
        }
    }
    
    std::cout << "Created " << particles.size() << " projectile particles\n";
    return particles;
}

// Function to create plate particles (rectangular block) - reduced resolution
std::vector<mysph::Particle<double>> create_plate() {
    std::vector<mysph::Particle<double>> particles;
    
    // Smaller plate dimensions
    const double plate_length = 0.005;    // Reduced from 0.025
    const double plate_width = 0.005;     // Reduced from 0.03
    const double plate_thickness = 0.0005; // 0.0025
    
    // Use a larger particle spacing for testing
    double test_dx = dx * 2.5;
    
    for (double x = 0; x <= plate_thickness; x += test_dx) {
        for (double y = 0; y <= plate_length; y += test_dx) {
            for (double z = -plate_width/2; z <= plate_width/2; z += test_dx) {
                mysph::Particle<double> p;
                p.r = {x, y, z};
                p.v = {0.0, 0.0, 0.0};
                p.m = test_dx * test_dx * test_dx * ro1;
                p.rho = p.rho0 = ro1;
                p.cs = cs1;
                p.material = 0; // plate
                p.G = G1;
                p.Yo = Yo1;
                p.k = E1;
                p.nu = 0.33; // Poisson ratio for aluminum
                particles.push_back(p);
            }
        }
    }
    
    std::cout << "Created " << particles.size() << " plate particles\n";
    return particles;
}

// Compute velocity gradient for a particle
void compute_velocity_gradient(std::vector<mysph::Particle<double>>& particles, int i, 
                               const std::vector<mysph::Particle<double>>& neighbors) {
    auto& pi = particles[i];
    pi.v00 = pi.v01 = pi.v02 = pi.v10 = pi.v11 = pi.v12 = pi.v20 = pi.v21 = pi.v22 = 0.0;
    
    for (auto& pj : neighbors) {
        auto rij = pi.r - pj.r;
        auto vij = pi.v - pj.v;
        auto grad_k = mysph::grad_kernel(rij, h);
        
        double factor = pj.m / pj.rho;
        
        pi.v00 += factor * vij[0] * grad_k[0];
        pi.v01 += factor * vij[0] * grad_k[1];
        pi.v02 += factor * vij[0] * grad_k[2];
        
        pi.v10 += factor * vij[1] * grad_k[0];
        pi.v11 += factor * vij[1] * grad_k[1];
        pi.v12 += factor * vij[1] * grad_k[2];
        
        pi.v20 += factor * vij[2] * grad_k[0];
        pi.v21 += factor * vij[2] * grad_k[1];
        pi.v22 += factor * vij[2] * grad_k[2];
    }
}

// Von Mises plasticity model
void apply_von_mises_plasticity(mysph::Particle<double>& p) {
    // Calculate stress deviator second invariant J2
    p.J2 = 0.5 * (p.s00*p.s00 + p.s11*p.s11 + p.s22*p.s22 + 
                   2*(p.s01*p.s01 + p.s02*p.s02 + p.s12*p.s12));
    
    // Check yield criterion
    if (p.J2 > (p.Yo * p.Yo / 3.0)) {
        // Calculate yield factor
        double yield_factor = p.Yo / (std::sqrt(3.0 * p.J2));
        
        // Scale down deviatoric stress components
        p.s00 *= yield_factor;
        p.s01 *= yield_factor;
        p.s02 *= yield_factor;
        p.s11 *= yield_factor;
        p.s12 *= yield_factor;
        p.s22 *= yield_factor;
    }
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

#include <fstream>
#include <iomanip>
#include <set>

// Function to write output in human-readable format
void write_evaluation_data(const std::string& filename, 
                          const std::vector<mysph::Particle<double>>& particles,
                          double time, int step) {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open output file: " << filename << std::endl;
        return;
    }
    
    // Write header
    file << "# SPH Impact Simulation - Time: " << time << " s, Step: " << step << "\n";
    file << "# Total particles: " << particles.size() << "\n\n";
    
    // Statistics section
    int plate_count = 0;
    int projectile_count = 0;
    double max_velocity = 0.0;
    double max_pressure = 0.0;
    double max_stress = 0.0;
    double max_plastic_strain = 0.0;
    double total_ke = 0.0;  // kinetic energy
    
    for (const auto& p : particles) {
        if (p.material == 0) plate_count++;
        else projectile_count++;
        
        double v_mag = mysph::abs(p.v);
        max_velocity = std::max(max_velocity, v_mag);
        max_pressure = std::max(max_pressure, std::abs(p.p));
        
        // von Mises equivalent stress
        double vm_stress = std::sqrt(0.5 * (
            std::pow(p.s00 - p.s11, 2) + 
            std::pow(p.s11 - p.s22, 2) + 
            std::pow(p.s22 - p.s00, 2) + 
            6.0 * (p.s01*p.s01 + p.s02*p.s02 + p.s12*p.s12)
        ));
        max_stress = std::max(max_stress, vm_stress);
        max_plastic_strain = std::max(max_plastic_strain, p.plastic_strain);
        
        // Kinetic energy
        total_ke += 0.5 * p.m * v_mag * v_mag;
    }
    
    file << "## Global Statistics\n";
    file << "Plate particles:        " << plate_count << "\n";
    file << "Projectile particles:   " << projectile_count << "\n";
    file << "Maximum velocity:       " << max_velocity << " m/s\n";
    file << "Maximum pressure:       " << max_pressure/1e6 << " MPa\n";
    file << "Maximum von Mises stress: " << max_stress/1e6 << " MPa\n";
    file << "Maximum plastic strain: " << max_plastic_strain << "\n";
    file << "Total kinetic energy:   " << total_ke << " J\n\n";
    
    // Sample particles section
    file << "## Sample Particle Data\n";
    file << std::setw(5) << "ID" << " | " 
          << std::setw(8) << "Material" << " | "
          << std::setw(10) << "X" << " | " 
          << std::setw(10) << "Y" << " | "
          << std::setw(10) << "Z" << " | "
          << std::setw(10) << "Vx" << " | "
          << std::setw(10) << "Vy" << " | "
          << std::setw(10) << "Vz" << " | "
          << std::setw(10) << "Density" << " | "
          << std::setw(10) << "Pressure" << " | "
          << std::setw(15) << "VM Stress\n";
    
    // Print header separator
    file << std::string(130, '-') << "\n";
    
    // Sample data (first 10 particles of each material and 10 particles with highest stress)
    int plate_samples = 0;
    int projectile_samples = 0;
    
    // Find high stress particles
    std::vector<std::pair<double, size_t>> stress_pairs;
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& p = particles[i];
        double vm_stress = std::sqrt(0.5 * (
            std::pow(p.s00 - p.s11, 2) + 
            std::pow(p.s11 - p.s22, 2) + 
            std::pow(p.s22 - p.s00, 2) + 
            6.0 * (p.s01*p.s01 + p.s02*p.s02 + p.s12*p.s12)
        ));
        stress_pairs.push_back({vm_stress, i});
    }
    
    // Sort by stress (descending)
    std::sort(stress_pairs.begin(), stress_pairs.end(), 
              [](const auto& a, const auto& b) { return a.first > b.first; });
    
    // Processed particles to avoid duplicates
    std::set<size_t> processed;
    
    // Print high stress particles (top 10)
    file << "# Top 10 particles by von Mises stress:\n";
    for (int i = 0; i < std::min(10, (int)stress_pairs.size()); i++) {
        size_t idx = stress_pairs[i].second;
        const auto& p = particles[idx];
        
        double vm_stress = std::sqrt(0.5 * (
            std::pow(p.s00 - p.s11, 2) + 
            std::pow(p.s11 - p.s22, 2) + 
            std::pow(p.s22 - p.s00, 2) + 
            6.0 * (p.s01*p.s01 + p.s02*p.s02 + p.s12*p.s12)
        ));
        
        file << std::setw(5) << idx << " | " 
             << std::setw(8) << (p.material == 0 ? "Plate" : "Projectile") << " | "
             << std::setw(10) << p.r[0] << " | " 
             << std::setw(10) << p.r[1] << " | "
             << std::setw(10) << p.r[2] << " | "
             << std::setw(10) << p.v[0] << " | "
             << std::setw(10) << p.v[1] << " | "
             << std::setw(10) << p.v[2] << " | "
             << std::setw(10) << p.rho << " | "
             << std::setw(10) << p.p/1e6 << " | "
             << std::setw(15) << vm_stress/1e6 << "\n";
             
        processed.insert(idx);
    }
    
    file << "\n# Sample particles from each material:\n";
    // Print sample particles from plate
    for (size_t i = 0; i < particles.size() && plate_samples < 10; i++) {
        const auto& p = particles[i];
        if (p.material == 0 && processed.find(i) == processed.end()) {
            double vm_stress = std::sqrt(0.5 * (
                std::pow(p.s00 - p.s11, 2) + 
                std::pow(p.s11 - p.s22, 2) + 
                std::pow(p.s22 - p.s00, 2) + 
                6.0 * (p.s01*p.s01 + p.s02*p.s02 + p.s12*p.s12)
            ));
            
            file << std::setw(5) << i << " | " 
                 << std::setw(8) << "Plate" << " | "
                 << std::setw(10) << p.r[0] << " | " 
                 << std::setw(10) << p.r[1] << " | "
                 << std::setw(10) << p.r[2] << " | "
                 << std::setw(10) << p.v[0] << " | "
                 << std::setw(10) << p.v[1] << " | "
                 << std::setw(10) << p.v[2] << " | "
                 << std::setw(10) << p.rho << " | "
                 << std::setw(10) << p.p/1e6 << " | "
                 << std::setw(15) << vm_stress/1e6 << "\n";
                 
            plate_samples++;
        }
    }
    
    // Print sample particles from projectile
    for (size_t i = 0; i < particles.size() && projectile_samples < 10; i++) {
        const auto& p = particles[i];
        if (p.material == 1 && processed.find(i) == processed.end()) {
            double vm_stress = std::sqrt(0.5 * (
                std::pow(p.s00 - p.s11, 2) + 
                std::pow(p.s11 - p.s22, 2) + 
                std::pow(p.s22 - p.s00, 2) + 
                6.0 * (p.s01*p.s01 + p.s02*p.s02 + p.s12*p.s12)
            ));
            
            file << std::setw(5) << i << " | " 
                 << std::setw(8) << "Projectile" << " | "
                 << std::setw(10) << p.r[0] << " | " 
                 << std::setw(10) << p.r[1] << " | "
                 << std::setw(10) << p.r[2] << " | "
                 << std::setw(10) << p.v[0] << " | "
                 << std::setw(10) << p.v[1] << " | "
                 << std::setw(10) << p.v[2] << " | "
                 << std::setw(10) << p.rho << " | "
                 << std::setw(10) << p.p/1e6 << " | "
                 << std::setw(15) << vm_stress/1e6 << "\n";
                 
            projectile_samples++;
        }
    }
    
    file.close();
}

// Function to write all particle properties to a file
void write_full_particle_data(const std::string& filename, 
                             const std::vector<mysph::Particle<double>>& particles,
                             double time, int step) {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open output file: " << filename << std::endl;
        return;
    }
    
    // Write header
    file << "# SPH Impact Simulation - Full Particle Data\n";
    file << "# Time: " << time << " s, Step: " << step << "\n";
    file << "# Total particles: " << particles.size() << "\n\n";
    
    // Column headers
    file << std::setw(5) << "ID" << " | " 
         << std::setw(8) << "Material" << " | "
         << std::setw(12) << "X" << " | " 
         << std::setw(12) << "Y" << " | "
         << std::setw(12) << "Z" << " | "
         << std::setw(12) << "Vx" << " | "
         << std::setw(12) << "Vy" << " | "
         << std::setw(12) << "Vz" << " | "
         << std::setw(12) << "Density" << " | "
         << std::setw(12) << "Pressure" << " | "
         << std::setw(12) << "Mass" << " | "
         << std::setw(12) << "S00" << " | "
         << std::setw(12) << "S01" << " | "
         << std::setw(12) << "S02" << " | "
         << std::setw(12) << "S11" << " | "
         << std::setw(12) << "S12" << " | "
         << std::setw(12) << "S22" << " | "
         << std::setw(12) << "VM_Stress" << " | "
         << std::setw(12) << "Pl_Strain" << "\n";
    
    // Print header separator
    file << std::string(250, '-') << "\n";
    
    // Print data for all particles
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& p = particles[i];
        
        // Calculate von Mises stress
        double vm_stress = std::sqrt(0.5 * (
            std::pow(p.s00 - p.s11, 2) + 
            std::pow(p.s11 - p.s22, 2) + 
            std::pow(p.s22 - p.s00, 2) + 
            6.0 * (p.s01*p.s01 + p.s02*p.s02 + p.s12*p.s12)
        ));
        
        file << std::scientific << std::setprecision(5);
        file << std::setw(5) << i << " | " 
             << std::setw(8) << (p.material == 0 ? "Plate" : "Projectile") << " | "
             << std::setw(12) << p.r[0] << " | " 
             << std::setw(12) << p.r[1] << " | "
             << std::setw(12) << p.r[2] << " | "
             << std::setw(12) << p.v[0] << " | "
             << std::setw(12) << p.v[1] << " | "
             << std::setw(12) << p.v[2] << " | "
             << std::setw(12) << p.rho << " | "
             << std::setw(12) << p.p << " | "
             << std::setw(12) << p.m << " | "
             << std::setw(12) << p.s00 << " | "
             << std::setw(12) << p.s01 << " | "
             << std::setw(12) << p.s02 << " | "
             << std::setw(12) << p.s11 << " | "
             << std::setw(12) << p.s12 << " | "
             << std::setw(12) << p.s22 << " | "
             << std::setw(12) << vm_stress << " | "
             << std::setw(12) << p.plastic_strain << "\n";
    }
    
    file.close();
}



int main() {
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
    
    // Time stepping
    int step = 0;
    double time = 0.0;
    
    std::filesystem::create_directory("output");
    
    while (time < max_time) {
        std::cout << "Step " << step << ", Time: " << time << " s\n";
        
        // Copy particles for integration
        auto next_particles = particles;
        
        // Density calculation (SPH summation)
        for (size_t i = 0; i < particles.size(); i++) {
            particles[i].rho = 0.0;
            
            for (size_t j = 0; j < particles.size(); j++) {
                particles[i].rho += particles[j].m * mysph::kernel(particles[i].r - particles[j].r, h);
            }

        }
        
        // Create neighbor lists for each particle
        std::vector<std::vector<mysph::Particle<double>>> neighbors(particles.size());
        for (size_t i = 0; i < particles.size(); i++) {
            for (size_t j = 0; j < particles.size(); j++) {
                if (i != j && mysph::abs(particles[i].r - particles[j].r) <= 2*h) {
                    neighbors[i].push_back(particles[j]);
                }
            }
        }
        
        // Compute velocity gradients
        for (size_t i = 0; i < particles.size(); i++) {
            compute_velocity_gradient(particles, i, neighbors[i]);
        }
        
        // Update pressure using EOS
        for (auto& p : particles) {
            compute_pressure(p);
        }
        
        // Update stress state
        for (auto& p : particles) {
            compute_stress_rate(p, dt);
        }
        
        // Compute acceleration (forces)
        for (size_t i = 0; i < particles.size(); i++) {
            auto& pi = particles[i];
            
            // Initialize force and viscous force
            pi.F = {0.0, 0.0, 0.0};
            pi.Fv = {0.0, 0.0, 0.0};
            
            for (auto& pj : neighbors[i]) {
                auto rij = pi.r - pj.r;
                auto vij = pi.v - pj.v;
                double r = mysph::abs(rij);
                
                if (r > 0) {
                    auto eij = rij * (1.0/r);
                    auto grad_k = mysph::grad_kernel(rij, h);
                    
                    // Artificial viscosity term
                    double dot_prod = vij[0]*rij[0] + vij[1]*rij[1] + vij[2]*rij[2];
                    double mu_ij = 0.0;
                    
                    if (dot_prod < 0) {
                        double c_ij = 0.5 * (pi.cs + pj.cs);
                        mu_ij = h * dot_prod / (r*r + avisc_eta*h*h);
                        double pi_ij = (-avisc_alpha * mu_ij * c_ij + avisc_beta * mu_ij * mu_ij) / 
                                    (0.5 * (pi.rho + pj.rho));
                        
                        // Viscous force contribution
                        pi.Fv = pi.Fv - grad_k * (pj.m * pi_ij);
                    }
                    
                    // Stress tensor force
                    double stress_term_i[3] = {
                        (pi.s00 + pi.p) * grad_k[0] + pi.s01 * grad_k[1] + pi.s02 * grad_k[2],
                        pi.s01 * grad_k[0] + (pi.s11 + pi.p) * grad_k[1] + pi.s12 * grad_k[2],
                        pi.s02 * grad_k[0] + pi.s12 * grad_k[1] + (pi.s22 + pi.p) * grad_k[2]
                    };
                    
                    double stress_term_j[3] = {
                        (pj.s00 + pj.p) * grad_k[0] + pj.s01 * grad_k[1] + pj.s02 * grad_k[2],
                        pj.s01 * grad_k[0] + (pj.s11 + pj.p) * grad_k[1] + pj.s12 * grad_k[2],
                        pj.s02 * grad_k[0] + pj.s12 * grad_k[1] + (pj.s22 + pj.p) * grad_k[2]
                    };
                    
                    // Total force contribution
                    double density_factor_i = pj.m / (pi.rho * pi.rho);
                    double density_factor_j = pj.m / (pj.rho * pj.rho);
                    
                    pi.F[0] -= density_factor_i * stress_term_i[0] + density_factor_j * stress_term_j[0];
                    pi.F[1] -= density_factor_i * stress_term_i[1] + density_factor_j * stress_term_j[1];
                    pi.F[2] -= density_factor_i * stress_term_i[2] + density_factor_j * stress_term_j[2];
                }
            }
        }
        
        // Integrate (update velocity and position)
        for (size_t i = 0; i < particles.size(); i++) {
            // V splitted (velocity update with viscous force)
            next_particles[i].vstar = particles[i].v + (particles[i].Fv) * (dt / particles[i].m);
            
            // Final velocity update with total force
            next_particles[i].v = next_particles[i].vstar + particles[i].F * (dt / particles[i].m);
            
            // Position update
            next_particles[i].r = particles[i].r + next_particles[i].v * dt;
        }
        
        // XSPH correction (optional)
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
        
        // Update particles
        particles = next_particles;

        // Write output
        if (step % 1 == 0) {
            // VTK file for visualization
            std::string vtk_filename = "output/impact-" + std::to_string(step) + ".vtp";
            write_particles_vtk(vtk_filename, particles);

            // Human-readable file for evaluation
            std::string eval_filename = "output/eval-" + std::to_string(step) + ".txt";
            write_full_particle_data(eval_filename, particles, time, step);
        }

        
        time += dt;
        step++;
    }
    
    return 0;
}
