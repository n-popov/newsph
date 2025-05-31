#pragma once

#include <fstream>
#include <iomanip>
#include <set>
#include <vector>
#include <cmath>
#include <algorithm> // For std::swap and std::abs
#include <iostream>  // For potential debugging
#include <iomanip>   // For formatting output

#include "particle.h"

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

// Small constant for floating point comparisons
const double JACOBI_EPSILON = 1e-10;
const int JACOBI_MAX_ROTATIONS = 50; // Max sweeps for 3x3, usually very few needed (e.g., 5-10)

// Computes eigenvalues and eigenvectors of a 3x3 real symmetric matrix.
std::pair<std::array<std::array<double, 3>, 3>, std::array<double, 3>> eigen_decomposition_3x3(const double S[3][3]) {
    std::array<std::array<double, 3>, 3> R_mat;
    std::array<double, 3> V;

    double A[3][3]; // Working copy of S
    // Initialize R_mat to identity matrix
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            A[i][j] = S[i][j];
            R_mat[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int sweep = 0; sweep < JACOBI_MAX_ROTATIONS; ++sweep) {
        // Calculate sum of squares of off-diagonal elements
        double sum_off_diag_sq = 0.0;
        sum_off_diag_sq += A[0][1] * A[0][1];
        sum_off_diag_sq += A[0][2] * A[0][2];
        sum_off_diag_sq += A[1][2] * A[1][2];
        sum_off_diag_sq *= 2.0; // Since A is symmetric

        if (sum_off_diag_sq < JACOBI_EPSILON * JACOBI_EPSILON) { // Arbitrary small threshold
            break; // Converged
        }

        // Perform rotations for (0,1), (0,2), (1,2)
        for (int p = 0; p < 3; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (std::abs(A[p][q]) < JACOBI_EPSILON / 100.0) { // Skip if element is already very small
                    continue;
                }

                double app = A[p][p];
                double aqq = A[q][q];
                double apq = A[p][q];

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
                // A[p][p] and A[q][q]
                A[p][p] = app - t * apq;
                A[q][q] = aqq + t * apq;
                A[p][q] = 0.0;
                A[q][p] = 0.0;

                // Other elements involving row/col p and q
                for (int k = 0; k < 3; ++k) {
                    if (k != p && k != q) {
                        double akp = A[k][p]; // Store old value before it's changed
                        double akq = A[k][q]; // Store old value
                        A[k][p] = c * akp - s * akq;
                        A[p][k] = A[k][p]; // Symmetry
                        A[k][q] = s * akp + c * akq;
                        A[q][k] = A[k][q]; // Symmetry
                    }
                }
                
                // Update Eigenvector matrix R_mat: R_new = R_old * G
                for (int k = 0; k < 3; ++k) {
                    double rkp = R_mat[k][p]; // Store old value
                    double rkq = R_mat[k][q]; // Store old value
                    R_mat[k][p] = c * rkp - s * rkq;
                    R_mat[k][q] = s * rkp + c * rkq;
                }
            }
        }
    }

    // Eigenvalues are the diagonal elements of the diagonalized A
    V[0] = A[0][0];
    V[1] = A[1][1];
    V[2] = A[2][2];

    // Sort eigenvalues and corresponding eigenvectors (simple bubble sort for 3 elements)
    // The eigenvectors in R_mat are columns.
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2 - i; ++j) {
            if (V[j] > V[j + 1]) {
                std::swap(V[j], V[j + 1]);
                // Swap corresponding columns in R_mat
                for (int k = 0; k < 3; ++k) {
                    std::swap(R_mat[k][j], R_mat[k][j + 1]);
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
std::array<std::array<double, 3>, 3> transform_diag_inv_3x3(const double rd[3], const std::array<std::array<double, 3>, 3> R_mat) {
    // Rab_ij = sum_k ( R_mat_ik * D_kk * (R_mat^T)_kj )
    // (R_mat^T)_kj = R_mat_jk
    // Rab_ij = sum_k ( R_mat_ik * rd_k * R_mat_jk )
    // Rab_ij = R_mat_i0 * rd_0 * R_mat_j0 + 
    //          R_mat_i1 * rd_1 * R_mat_j1 +
    //          R_mat_i2 * rd_2 * R_mat_j2

    std::array<std::array<double, 3>, 3> Rab;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Rab[i][j] = R_mat[i][0] * rd[0] * R_mat[j][0] +
                        R_mat[i][1] * rd[1] * R_mat[j][1] +
                        R_mat[i][2] * rd[2] * R_mat[j][2];
        }
    }

    return Rab;
}