#pragma once

#include <cmath>
#include <numeric>
#include <numbers>
#include <functional>
#include <concepts>
#include <stdexcept>

#include "../utils/common.h"

enum class Kernel {
    TUTORIAL,
    GAUSSIAN,
    CUBIC_SPLINE,
    WENDLAND_QUINTIC
};

namespace mysph {

namespace {
    // --- Tutorial Kernel ---
    template <class T>
    T qsmall(T q) {
        return T(6) * (std::pow(q, 3) - std::pow(q, 2)) + T(1);
    }

    template <class T>
    T qbig(T q) {
        return T(2) * std::pow((T(1) - q), 3);
    }

    inline double sigma3(double h) {
        return 8 / (std::numbers::pi * h * h * h);
    }

    template <class T>
    double _kernel(T r, double h) {
        auto q = abs(r) / h;
        if (2 * q <= 1) {
            return qsmall(q);
        }
        if (2 * q <= 2) {
            return qbig(q);
        }
        return 0;
    }

    template <class T>
    double _grad_kernel(T r, double h) {
        auto q = abs(r) / h;
        if (2 * q <= 1) {
            return 18 * q * q - 12 * q;
        }
        if (2 * q <= 2) {
            return -6 * (q - 1) * (q - 1);
        }
        return 0;
    }

    template <class T>
    double kernel_tutorial(T r, double h) {
        return _kernel(r, h) * sigma3(h);
    }

    template <class T>
    T grad_kernel_tutorial(T r, double h) {
        auto den = h * abs(r);
        if (den < 1e-12) return T();
        return r * (_grad_kernel(r, h) * (sigma3(h)) / den);
    }

    // --- Gaussian Kernel ---
    double sigma_gaussian(double h) {
        return 1 / (std::pow(std::numbers::pi, 1.5) * h * h * h);
    }

    double gaussian(vec3<double> r, double h) {
        auto q = abs(r) / h;
        if (q < 3) {
            return sigma_gaussian(h) * std::exp(-q * q);
        } else {
            return 0.;
        }
    }

    vec3<double> grad_gaussian(vec3<double> r, double h) {
        auto r_norm = abs(r);
        auto den = h * r_norm;
        if (den < 1e-12) return {};

        auto q = r_norm / h;
        if (q < 3) {
            return r * (sigma_gaussian(h) * (-2) * std::exp(-q * q) / (h * h));
        } else {
            return {};
        }
    }

    // --- Cubic Spline Kernel (3D) ---
    double sigma_cubic_spline(double h) {
        return 1.0 / (std::numbers::pi * h * h * h);
    }

    double cubic_spline(vec3<double> r, double h) {
        const auto q = abs(r) / h;
        const auto sigma = sigma_cubic_spline(h);

        if (q <= 1.0) {
            return sigma * (1.0 - 1.5 * q * q * (1.0 - 0.5 * q));
        } else if (q <= 2.0) {
            const auto term = 2.0 - q;
            return sigma * 0.25 * term * term * term;
        }
        return 0.0;
    }

    vec3<double> grad_cubic_spline(vec3<double> r, double h) {
        const auto r_norm = abs(r);
        const auto den = h * r_norm;
        if (den < 1e-12) return {};

        const auto q = r_norm / h;
        const auto sigma = sigma_cubic_spline(h);
        double factor = 0.0;

        if (q <= 1.0) {
            // d/dq (1 - 1.5q^2 + 0.75q^3) = -3q + 2.25q^2
            factor = -3.0 * q + 2.25 * q * q;
        } else if (q <= 2.0) {
            // d/dq (0.25 * (2-q)^3) = -0.75 * (2-q)^2
            const auto term = 2.0 - q;
            factor = -0.75 * term * term;
        }

        return r * (sigma * factor / den);
    }

    // --- Wendland Quintic Kernel (C2, 3D) ---
    double alpha_wendland_quintic(double h) {
        return 21.0 / (16.0 * std::numbers::pi * h * h * h);
    }

    double wendland_quintic(vec3<double> r, double h) {
        const auto q = abs(r) / h;

        if (q <= 2.0) {
            const auto alpha = alpha_wendland_quintic(h);
            const auto term1 = 1.0 - 0.5 * q;
            const auto term2 = 2.0 * q + 1.0;
            return alpha * (term1 * term1 * term1 * term1) * term2;
        }
        return 0.0;
    }

    vec3<double> grad_wendland_quintic(vec3<double> r, double h) {
        const auto r_norm = abs(r);
        const auto den = h * r_norm;
        if (den < 1e-12) return {};

        const auto q = r_norm / h;
        if (q <= 2.0) {
            const auto alpha = alpha_wendland_quintic(h);
            const auto term = 1.0 - 0.5 * q;
            // dW/dq = alpha * -5q * (1-q/2)^3
            const double factor = -5.0 * q * (term * term * term);
            return r * (alpha * factor / den);
        }
        return {};
    }

}; // end anonymous namespace

double kernel(vec3<double> r, double h, Kernel k_type) {
    if (k_type == Kernel::TUTORIAL)         return kernel_tutorial(r, h);
    if (k_type == Kernel::GAUSSIAN)         return gaussian(r, h);
    if (k_type == Kernel::CUBIC_SPLINE)     return cubic_spline(r, h);
    if (k_type == Kernel::WENDLAND_QUINTIC) return wendland_quintic(r, h);

    throw std::runtime_error("unspecified kernel");
}

vec3<double> grad_kernel(vec3<double> r, double h, Kernel k_type) {
    if (k_type == Kernel::TUTORIAL)         return grad_kernel_tutorial(r, h);
    if (k_type == Kernel::GAUSSIAN)         return grad_gaussian(r, h);
    if (k_type == Kernel::CUBIC_SPLINE)     return grad_cubic_spline(r, h);
    if (k_type == Kernel::WENDLAND_QUINTIC) return grad_wendland_quintic(r, h);

    throw std::runtime_error("unspecified kernel");
}

}; // end namespace mysph
