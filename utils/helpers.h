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
