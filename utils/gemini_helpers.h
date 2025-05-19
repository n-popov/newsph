#pragma once

#include <vector>
#include <cmath>
#include <algorithm> // For std::swap and std::abs
#include <iostream>  // For potential debugging
#include <iomanip>   // For formatting output

#include "physics.h"

// Small constant for floating point comparisons
const double JACOBI_EPSILON = 1e-10;
const int JACOBI_MAX_ROTATIONS = 50; // Max sweeps for 3x3, usually very few needed (e.g., 5-10)

// Helper function to print a 3x3 matrix (for debugging)
void print_matrix3x3(const double M[3][3], const std::string& name) {
    std::cout << name << " = [\n";
    for (int i = 0; i < 3; ++i) {
        std::cout << "  ";
        for (int j = 0; j < 3; ++j) {
            std::cout << std::fixed << std::setw(12) << std::setprecision(6) << M[i][j] << (j == 2 ? "" : ", ");
        }
        std::cout << (i == 2 ? "\n]" : ";\n");
    }
    std::cout << std::endl;
}

// Helper function to print a 3-vector (for debugging)
void print_vector3(const double V[3], const std::string& name) {
    std::cout << name << " = [";
    for (int i = 0; i < 3; ++i) {
        std::cout << std::fixed << std::setw(12) << std::setprecision(6) << V[i] << (i == 2 ? "" : ", ");
    }
    std::cout << "]\n" << std::endl;
}


/**
 * @brief Computes eigenvalues and eigenvectors of a 3x3 real symmetric matrix.
 *
 * Uses a Jacobi-like method to diagonalize the matrix.
 * Eigenvalues are sorted in ascending order.
 * Eigenvectors are stored as columns in R_mat.
 *
 * @param S Input 3x3 symmetric matrix.
 * @param R_mat Output 3x3 matrix whose columns are the eigenvectors.
 * @param V Output array of 3 eigenvalues, sorted in ascending order.
 */
void eigen_decomposition_3x3(const double S[3][3], double R_mat[3][3], double V[3]) {
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
}


/**
 * @brief Transforms a diagonal matrix (represented by rd) back to the original frame.
 *
 * Computes Rab = R_mat * D * R_mat^T, where D is a diagonal matrix
 * with elements from rd. R_mat contains orthonormal eigenvectors as columns.
 *
 * @param rd Input array of 3 diagonal elements.
 * @param R_mat Input 3x3 rotation matrix (eigenvectors as columns).
 * @param Rab Output 3x3 transformed symmetric matrix.
 */
void transform_diag_inv_3x3(const double rd[3], const double R_mat[3][3], double Rab[3][3]) {
    // Rab_ij = sum_k ( R_mat_ik * D_kk * (R_mat^T)_kj )
    // (R_mat^T)_kj = R_mat_jk
    // Rab_ij = sum_k ( R_mat_ik * rd_k * R_mat_jk )
    // Rab_ij = R_mat_i0 * rd_0 * R_mat_j0 + 
    //          R_mat_i1 * rd_1 * R_mat_j1 +
    //          R_mat_i2 * rd_2 * R_mat_j2

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Rab[i][j] = R_mat[i][0] * rd[0] * R_mat[j][0] +
                        R_mat[i][1] * rd[1] * R_mat[j][1] +
                        R_mat[i][2] * rd[2] * R_mat[j][2];
        }
    }
}


// --- Example Usage (Place in a main function or test suite) ---
/*
int main() {
    // Example symmetric matrix S
    double S[3][3] = {
        {5.0, 2.0, 0.0},
        {2.0, 6.0, -1.0},
        {0.0, -1.0, 7.0}
    };
    // S should be: (using Octave/Matlab `eig([5 2 0; 2 6 -1; 0 -1 7])`)
    // Eigenvalues approx: 3.649, 6.000, 8.351
    // Eigenvectors are columns of R_mat.

    double R_mat[3][3];
    double V[3];
    double Rab[3][3];

    std::cout << "Original Matrix S:" << std::endl;
    print_matrix3x3(S, "S");

    eigen_decomposition_3x3(S, R_mat, V);

    std::cout << "Eigenvalues (V):" << std::endl;
    print_vector3(V, "V");
    std::cout << "Eigenvectors (columns of R_mat):" << std::endl;
    print_matrix3x3(R_mat, "R_mat");

    // Test: Reconstruct S from R_mat, V. Should be S = R_mat * diag(V) * R_mat^T
    // This is what transform_diag_inv_3x3 does.
    transform_diag_inv_3x3(V, R_mat, Rab);
    std::cout << "Reconstructed S (R_mat * diag(V) * R_mat^T):" << std::endl;
    print_matrix3x3(Rab, "Rab_reconstructed");

    // Example for transform_diag_inv_3x3 with some artificial rd
    double rd_example[3] = {1.0, -0.5, 2.0}; // Artificial diagonal elements
    transform_diag_inv_3x3(rd_example, R_mat, Rab); // R_mat comes from S's decomposition
    
    std::cout << "Transformed Rab (using rd_example and R_mat from S):" << std::endl;
    print_matrix3x3(Rab, "Rab_example");

    return 0;
}
*/

// Forward declaration if needed, or include the Eigen library if you implement it yourself
// For Monaghan Artificial Stress
// void eigen_decomposition_3x3(const double S[3][3], double R_mat[3][3], double V[3]);
// void transform_diag_inv_3x3(const double rd[3], const double R_mat[3][3], double Rab[3][3]);

// Function to compute Monaghan Artificial Stress
void compute_monaghan_artificial_stress(mysph::Particle<double>& p, double eps) {
    double rhoi = p.rho;
    if (rhoi < 1e-9) { // Avoid division by zero if density is tiny
        p.as00 = 0.0; p.as01 = 0.0; p.as02 = 0.0;
        p.as11 = 0.0; p.as12 = 0.0; p.as22 = 0.0;
        return;
    }
    double rhoi21 = 1.0 / (rhoi * rhoi);

    double S[3][3]; // Total stress tensor
    double R_mat[3][3]; // Matrix of Eigenvectors (columns)
    double V[3];      // Eigenvalues (principal stresses)
    double rd[3];     // Artificial stress components in principal directions
    double Rab[3][3]; // Artificial stress tensor in original coordinates

    // Construct total stress tensor S = s_dev - p * I
    // Assuming p.pressure is positive for compression
    S[0][0] = p.s00 - p.p; S[0][1] = p.s01;        S[0][2] = p.s02;
    S[1][0] = p.s01;        S[1][1] = p.s11 - p.p; S[1][2] = p.s12;
    S[2][0] = p.s02;        S[2][1] = p.s12;        S[2][2] = p.s22 - p.p;

    // Compute the principal stresses and eigenvectors
    // You'll need a 3x3 symmetric eigenvalue solver here.
    // For simplicity, let's assume you have:
    eigen_decomposition_3x3(S, R_mat, V); // Fills R_mat and V

    // Artificial stress corrections in principal directions
    rd[0] = (V[0] > 0) ? -eps * V[0] * rhoi21 : 0.0;
    rd[1] = (V[1] > 0) ? -eps * V[1] * rhoi21 : 0.0;
    rd[2] = (V[2] > 0) ? -eps * V[2] * rhoi21 : 0.0;

    // Transform artificial stresses back to original frame: Rab = R_mat * diag(rd) * R_mat_transpose
    transform_diag_inv_3x3(rd, R_mat, Rab); // Fills Rab

    // Store the values in the particle (assuming p has r00, r01 etc. members)
    p.as00 = Rab[0][0];
    p.as01 = Rab[0][1]; // or Rab[1][0] due to symmetry
    p.as02 = Rab[0][2]; // or Rab[2][0]
    p.as11 = Rab[1][1];
    p.as12 = Rab[1][2]; // or Rab[2][1]
    p.as22 = Rab[2][2];
}


// Modified compute_stress_rate
void compute_stress_rate_and_artificial_terms(mysph::Particle<double>& p, double dt, double artificial_stress_eps) {
    // Calculate deviatoric strain rate
    double div_v = p.v00 + p.v11 + p.v22;
    double eps00 = p.v00 - div_v/3.0;
    double eps01 = 0.5 * (p.v01 + p.v10);
    double eps02 = 0.5 * (p.v02 + p.v20);
    double eps11 = p.v11 - div_v/3.0;
    double eps12 = 0.5 * (p.v12 + p.v21);
    double eps22 = p.v22 - div_v/3.0;

    // Update deviatoric stress using Hooke's law (elastic trial stress)
    p.s00 += 2 * p.G * eps00 * dt;
    p.s01 += 2 * p.G * eps01 * dt;
    p.s02 += 2 * p.G * eps02 * dt;
    p.s11 += 2 * p.G * eps11 * dt;
    p.s12 += 2 * p.G * eps12 * dt;
    p.s22 += 2 * p.G * eps22 * dt;

    // Apply plasticity (modifies p.s00, p.s01, etc.)
    apply_von_mises_plasticity(p);

    // -------- NEWLY ADDED --------
    // Compute Monaghan Artificial Stress based on the *updated* deviatoric stresses and current pressure
    // This requires p.pressure and p.rho to be up-to-date.
    if (artificial_stress_eps > 0) { // Only compute if eps is non-zero
        compute_monaghan_artificial_stress(p, artificial_stress_eps);
    } else {
        // Ensure artificial stress terms are zero if not computed
        p.as00 = 0.0; p.as01 = 0.0; p.as02 = 0.0;
        p.as11 = 0.0; p.as12 = 0.0; p.as22 = 0.0;
    }
    // -------- END NEWLY ADDED --------
}
