#pragma once

#include <cmath>
#include <numeric>
#include <numbers>
#include <functional>
#include <concepts>

#include "../utils/common.h"

namespace mysph {

namespace {
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
};

template <class T>
double kernel(T r, double h) {
    return _kernel(r, h) * sigma3(h);
}

template <class T>
T grad_kernel(T r, double h) {
    auto den = h * abs(r);

    if (den == 0) return T();
    
    return r * (_grad_kernel(r, h) * (sigma3(h)) / den);
}

};