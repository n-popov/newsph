#pragma once

#include <array>
#include <numeric>
#include <limits>
#include <fstream>
#include <iomanip>
#include <set>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <future>
#include <thread>
#include <algorithm>

#include "../config/simulation_config.h"
#include "../method/particle.h"

namespace {
    struct GroupStats {
            size_t count = 0;
            double total_mass = 0.0;
            
            // Density statistics
            double rho_min = std::numeric_limits<double>::max();
            double rho_max = std::numeric_limits<double>::lowest();
            double rho_avg = 0.0;
            
            // Velocity statistics
            double vel_mag_min = std::numeric_limits<double>::max();
            double vel_mag_max = std::numeric_limits<double>::lowest();
            double vel_mag_median = 0.0;
            double vel_x_avg = 0.0;
            double vel_y_avg = 0.0;
            double vel_z_avg = 0.0;
            
            // Position statistics
            double x_min = std::numeric_limits<double>::max();
            double x_max = std::numeric_limits<double>::lowest();
            double y_min = std::numeric_limits<double>::max();
            double y_max = std::numeric_limits<double>::lowest();
            double z_min = std::numeric_limits<double>::max();
            double z_max = std::numeric_limits<double>::lowest();
            double com_x = 0.0; // Center of mass
            double com_y = 0.0;
            double com_z = 0.0;
            
            // Force statistics
            double F_mag_min = std::numeric_limits<double>::max();
            double F_mag_max = std::numeric_limits<double>::lowest();
            double Fx_min = std::numeric_limits<double>::max();
            double Fx_max = std::numeric_limits<double>::lowest();
            double Fy_min = std::numeric_limits<double>::max();
            double Fy_max = std::numeric_limits<double>::lowest();
            double Fz_min = std::numeric_limits<double>::max();
            double Fz_max = std::numeric_limits<double>::lowest();
            
            // Energy statistics
            double total_kinetic = 0.0;
            double total_internal = 0.0;
            
            // Stress statistics
            double J2_min = std::numeric_limits<double>::max();
            double J2_max = std::numeric_limits<double>::lowest();
            double J2_avg = 0.0;
        };
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
         << std::setw(12) << "Pl_Strain" << " | "
         << std::setw(12) << "Fv_x" << " | "
         << std::setw(12) << "Fv_y" << " | "
         << std::setw(12) << "Fv_z" << " | "
         << std::setw(12) << "F_x" << " | "
         << std::setw(12) << "F_y" << " | "
         << std::setw(12) << "F_z" << "\n";
    
    // Print header separator
    file << std::string(350, '-') << "\n";
    
    // Print data for all particles
    for (size_t i = 0; i < particles.size(); i++) {
        const auto& p = particles[i];
        
        // Calculate von Mises stress
        double vm_stress = std::sqrt(0.5 * (
            std::pow(p.stress[0] - p.stress[4], 2) + 
            std::pow(p.stress[4] - p.stress[8], 2) + 
            std::pow(p.stress[8] - p.stress[0], 2) + 
            6.0 * (p.stress[1]*p.stress[1] + p.stress[2]*p.stress[2] + p.stress[5]*p.stress[7])
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
             << std::setw(12) << p.stress[0] << " | "
             << std::setw(12) << p.stress[1] << " | "
             << std::setw(12) << p.stress[2] << " | "
             << std::setw(12) << p.stress[4] << " | "
             << std::setw(12) << p.stress[5] << " | "
             << std::setw(12) << p.stress[8] << " | "
             << std::setw(12) << vm_stress << " | "
             << std::setw(12) << p.plastic_strain << " | "
             << std::setw(12) << p.Fv[0] << " | "
             << std::setw(12) << p.Fv[1] << " | "
             << std::setw(12) << p.Fv[2] << " | "
             << std::setw(12) << p.F[0] << " | "
             << std::setw(12) << p.F[1] << " | "
             << std::setw(12) << p.F[2] << "\n";
    }
    
    file.close();
}

void write_integral_characteristics(const std::string& filename, 
                             const std::vector<mysph::Particle<double>>& particles,
                             double time, int step) {
    std::ofstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Failed to open output file: " << filename << std::endl;
        return;
    }
    
    // Write header
    file << "# SPH Impact Simulation - Integral Characteristics\n";
    file << "# Time: " << time << " s, Step: " << step << "\n";
    file << "# Total particles: " << particles.size() << "\n\n";

    // Lambda to compute group statistics
    auto compute_stats = [](const std::vector<mysph::Particle<double>>& particles, int material) {
        
        
        GroupStats stats;
        std::vector<double> vel_magnitudes;
        std::vector<double> J2_values;
        
        for (const auto& p : particles) {
            if (p.is_fake || p.material != material) continue;
            
            stats.count++;
            stats.total_mass += p.m;
            
            // Velocity magnitude
            const double v_mag = std::sqrt(p.v[0]*p.v[0] + p.v[1]*p.v[1] + p.v[2]*p.v[2]);
            vel_magnitudes.push_back(v_mag);
            
            // Density stats
            stats.rho_min = std::min(stats.rho_min, p.rho);
            stats.rho_max = std::max(stats.rho_max, p.rho);
            stats.rho_avg += p.rho;
            
            // Position stats
            stats.x_min = std::min(stats.x_min, p.r[0]);
            stats.x_max = std::max(stats.x_max, p.r[0]);
            stats.y_min = std::min(stats.y_min, p.r[1]);
            stats.y_max = std::max(stats.y_max, p.r[1]);
            stats.z_min = std::min(stats.z_min, p.r[2]);
            stats.z_max = std::max(stats.z_max, p.r[2]);
            
            // Center of mass
            stats.com_x += p.m * p.r[0];
            stats.com_y += p.m * p.r[1];
            stats.com_z += p.m * p.r[2];
            
            // Force stats
            const double F_mag = std::sqrt(p.F[0]*p.F[0] + p.F[1]*p.F[1] + p.F[2]*p.F[2]);
            stats.F_mag_min = std::min(stats.F_mag_min, F_mag);
            stats.F_mag_max = std::max(stats.F_mag_max, F_mag);
            stats.Fx_min = std::min(stats.Fx_min, p.F[0]);
            stats.Fx_max = std::max(stats.Fx_max, p.F[0]);
            stats.Fy_min = std::min(stats.Fy_min, p.F[1]);
            stats.Fy_max = std::max(stats.Fy_max, p.F[1]);
            stats.Fz_min = std::min(stats.Fz_min, p.F[2]);
            stats.Fz_max = std::max(stats.Fz_max, p.F[2]);
            
            // Velocity components for average
            stats.vel_x_avg += p.v[0];
            stats.vel_y_avg += p.v[1];
            stats.vel_z_avg += p.v[2];
            
            // Energy calculations
            stats.total_kinetic += 0.5 * p.m * (p.v[0]*p.v[0] + p.v[1]*p.v[1] + p.v[2]*p.v[2]);
            stats.total_internal += p.m * p.e;
            
            // Stress invariant J2
            stats.J2_min = std::min(stats.J2_min, p.J2);
            stats.J2_max = std::max(stats.J2_max, p.J2);
            stats.J2_avg += p.J2;
            J2_values.push_back(p.J2);
        }
        
        // Post-process calculations
        if (stats.count > 0) {
            // Averages
            stats.rho_avg /= stats.count;
            stats.vel_x_avg /= stats.count;
            stats.vel_y_avg /= stats.count;
            stats.vel_z_avg /= stats.count;
            stats.J2_avg /= stats.count;
            
            // Center of mass
            stats.com_x /= stats.total_mass;
            stats.com_y /= stats.total_mass;
            stats.com_z /= stats.total_mass;
            
            // Velocity magnitude stats
            if (!vel_magnitudes.empty()) {
                std::sort(vel_magnitudes.begin(), vel_magnitudes.end());
                stats.vel_mag_min = vel_magnitudes.front();
                stats.vel_mag_max = vel_magnitudes.back();
                stats.vel_mag_median = vel_magnitudes.size() % 2 == 0 ?
                    (vel_magnitudes[vel_magnitudes.size()/2 - 1] + vel_magnitudes[vel_magnitudes.size()/2]) * 0.5 :
                    vel_magnitudes[vel_magnitudes.size()/2];
            }
        }
        
        return stats;
    };
    
    // Compute statistics for plate (material=0) and projectile (material=1)
    const auto plate_stats = compute_stats(particles, 0);
    const auto proj_stats = compute_stats(particles, 1);
    
    // Helper lambda for writing group statistics
    auto write_group = [&](const std::string& name, const GroupStats& stats) {
        file << "[" << name << "]\n";
        file << "count = " << stats.count << "\n";
        file << "total_mass = " << stats.total_mass << "\n";
        file << "density_min = " << stats.rho_min << "\n";
        file << "density_max = " << stats.rho_max << "\n";
        file << "density_avg = " << stats.rho_avg << "\n";
        file << "velocity_mag_min = " << stats.vel_mag_min << "\n";
        file << "velocity_mag_max = " << stats.vel_mag_max << "\n";
        file << "velocity_mag_median = " << stats.vel_mag_median << "\n";
        file << "velocity_x_avg = " << stats.vel_x_avg << "\n";
        file << "velocity_y_avg = " << stats.vel_y_avg << "\n";
        file << "velocity_z_avg = " << stats.vel_z_avg << "\n";
        file << "x_min = " << stats.x_min << "\n";
        file << "x_max = " << stats.x_max << "\n";
        file << "y_min = " << stats.y_min << "\n";
        file << "y_max = " << stats.y_max << "\n";
        file << "z_min = " << stats.z_min << "\n";
        file << "z_max = " << stats.z_max << "\n";
        file << "com_x = " << stats.com_x << "\n";
        file << "com_y = " << stats.com_y << "\n";
        file << "com_z = " << stats.com_z << "\n";
        file << "force_mag_min = " << stats.F_mag_min << "\n";
        file << "force_mag_max = " << stats.F_mag_max << "\n";
        file << "force_x_min = " << stats.Fx_min << "\n";
        file << "force_x_max = " << stats.Fx_max << "\n";
        file << "force_y_min = " << stats.Fy_min << "\n";
        file << "force_y_max = " << stats.Fy_max << "\n";
        file << "force_z_min = " << stats.Fz_min << "\n";
        file << "force_z_max = " << stats.Fz_max << "\n";
        file << "kinetic_energy = " << stats.total_kinetic << "\n";
        file << "internal_energy = " << stats.total_internal << "\n";
        file << "J2_stress_min = " << stats.J2_min << "\n";
        file << "J2_stress_max = " << stats.J2_max << "\n";
        file << "J2_stress_avg = " << stats.J2_avg << "\n\n";
    };
    
    // Write results for both groups
    write_group("Plate", plate_stats);
    write_group("Projectile", proj_stats);
}


// Small constant for floating point comparisons
const double JACOBI_EPSILON = 1e-10;
const int JACOBI_MAX_ROTATIONS = 50; // Max sweeps for 3x3, usually very few needed (e.g., 5-10)

// Computes eigenvalues and eigenvectors of a 3x3 real symmetric matrix.
std::pair<std::array<double, 9>, std::array<double, 3>> 
eigen_decomposition_3x3(const std::array<double, 9>& S) {
    std::array<double, 9> R_mat; // Flattened 3x3 eigenvector matrix
    std::array<double, 3> V;     // Eigenvalues
    std::array<double, 9> A;     // Working copy of S

    // Initialize R_mat to identity matrix and copy S to A
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A[i*3 + j] = S[i*3 + j];
            R_mat[i*3 + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int sweep = 0; sweep < JACOBI_MAX_ROTATIONS; ++sweep) {
        // Calculate sum of squares of off-diagonal elements
        double sum_off_diag_sq = 0.0;
        sum_off_diag_sq += A[1] * A[1]; // A[0][1]
        sum_off_diag_sq += A[2] * A[2]; // A[0][2]
        sum_off_diag_sq += A[5] * A[5]; // A[1][2]
        sum_off_diag_sq *= 2.0; // Since A is symmetric

        if (sum_off_diag_sq < JACOBI_EPSILON * JACOBI_EPSILON) {
            break; // Converged
        }

        // Perform rotations for (0,1), (0,2), (1,2)
        for (int p = 0; p < 3; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                double apq = A[p*3 + q];
                if (std::abs(apq) < JACOBI_EPSILON / 100.0) {
                    continue;
                }

                double app = A[p*3 + p];
                double aqq = A[q*3 + q];
                double tau = (aqq - app) / (2.0 * apq);
                double t;
                if (tau >= 0.0) {
                    t = 1.0 / (tau + std::sqrt(1.0 + tau * tau));
                } else {
                    t = -1.0 / (-tau + std::sqrt(1.0 + tau * tau));
                }
                double c = 1.0 / std::sqrt(1.0 + t * t);
                double s = t * c;

                // Update A matrix: A_new = G^T * A * G
                A[p*3 + p] = app - t * apq;
                A[q*3 + q] = aqq + t * apq;
                A[p*3 + q] = 0.0;
                A[q*3 + p] = 0.0;

                // Other elements involving row/col p and q
                for (int k = 0; k < 3; ++k) {
                    if (k != p && k != q) {
                        double akp = A[k*3 + p]; // Store old value
                        double akq = A[k*3 + q]; // Store old value
                        A[k*3 + p] = c * akp - s * akq;
                        A[p*3 + k] = A[k*3 + p]; // Symmetry
                        A[k*3 + q] = s * akp + c * akq;
                        A[q*3 + k] = A[k*3 + q]; // Symmetry
                    }
                }
                
                // Update Eigenvector matrix R_mat: R_new = R_old * G
                for (int k = 0; k < 3; ++k) {
                    double rkp = R_mat[k*3 + p]; // Store old value
                    double rkq = R_mat[k*3 + q]; // Store old value
                    R_mat[k*3 + p] = c * rkp - s * rkq;
                    R_mat[k*3 + q] = s * rkp + c * rkq;
                }
            }
        }
    }

    V[0] = A[0];
    V[1] = A[4];
    V[2] = A[8];

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2 - i; ++j) {
            if (V[j] > V[j + 1]) {
                std::swap(V[j], V[j + 1]);
                for (int k = 0; k < 3; ++k) {
                    std::swap(R_mat[k*3 + j], R_mat[k*3 + j + 1]);
                }
            }
        }
    }

    return {R_mat, V};
}


/**
 * Computes Rab = R_mat * D * R_mat^T, where D is a diagonal matrix
 * with elements from rd. R_mat contains orthonormal eigenvectors as columns.
 */
std::array<double, 9> transform_diag_inv_3x3(const std::array<double, 3>& rd, 
                                            const std::array<double, 9>& R_mat) {
    // Computes R * diag(rd) * R^T 
    // Rab_ij = sum_k ( R_mat_ik * rd_k * R_mat_jk )
    
    std::array<double, 9> Rab;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Rab[i*3 + j] = R_mat[i*3 + 0] * rd[0] * R_mat[j*3 + 0] +
                           R_mat[i*3 + 1] * rd[1] * R_mat[j*3 + 1] +
                           R_mat[i*3 + 2] * rd[2] * R_mat[j*3 + 2];
        }
    }

    return Rab;
}

const config::MaterialProperties& get_material_properties(const mysph::Particle<double>& particle, const config::SimulationConfig& config) {
    if (particle.material == 0) {
        return config.aluminum_props;
    } else if (particle.material == 1) {
        return config.steel_props;
    }

    throw std::runtime_error("unknown material of particle");
}

template<typename Function, typename... Args>
void parallelize(bool is_enabled,
                 Function&& func, 
                 std::vector<mysph::Particle<double>>& particles, 
                 std::vector<std::vector<mysph::Particle<double>*>>& neighbors, 
                 Args&&... args) {
    if (!is_enabled) {
        for (auto i = 0u; i < std::size(particles); i++) {
            func(particles[i], neighbors[i], args...);
        }

        return;
    }

    const auto num_threads = 4u;
    const auto num_particles = particles.size();
    const auto chunk_size = (num_particles + num_threads - 1) / num_threads;

    std::vector<std::future<void>> futures;

    for (auto thread_id = 0u; thread_id < num_threads; ++thread_id) {
        const auto start = thread_id * chunk_size;
        const auto end = std::min(start + chunk_size, num_particles);

        if (start >= end) {
            continue;
        }

        futures.emplace_back(std::async(std::launch::async, 
            [&func, &particles, &neighbors, start, end, &args...]() {
                for (size_t i = start; i < end; ++i) {
                    func(particles[i], neighbors[i], args...);
                }
            }));
    }

    for (auto& future : futures) {
        future.get();
    }
}