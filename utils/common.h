#pragma once

#include <array>
#include <cmath>
#include <numeric>

namespace mysph {

constexpr const std::array<double, 9> SINGILAR9 = {1., 0., 0., 0., 1., 0., 0., 0., 1.};

using vector3d = std::array<double, 3>;

template<typename T, size_t N>
std::array<T, N> operator+(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    std::array<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] + rhs[i];
    }
    return result;
}

template<typename T, size_t N>
std::array<T, N> operator-(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    std::array<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] - rhs[i];
    }
    return result;
}

template<typename T, size_t N>
std::array<T, N> operator*(const std::array<T, N>& lhs, T rhs) {
    std::array<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] * rhs;
    }
    return result;
}

template<typename T, size_t N>
std::array<T, N> operator/(const std::array<T, N>& lhs, T rhs) {
    std::array<T, N> result;
    for (size_t i = 0; i < N; ++i) {
        result[i] = lhs[i] / rhs;
    }
    return result;
}

template<typename T, size_t N>
T abs(const std::array<T, N>& v) {
    T sum = 0;
    for (const auto& val : v) {
        sum += val * val;
    }
    return std::sqrt(sum);
}

template<typename T, size_t N>
std::array<T, N> normalize(const std::array<T, N>& v) {
    T length = abs(v);
    if (length == 0) return v;
    return v / length;
}

template<typename T, size_t N>
T operator*(const std::array<T, N>& lhs, const std::array<T, N>& rhs) {
    return std::inner_product(std::begin(lhs), std::end(lhs), std::begin(rhs), T(0));
}

template<typename T>
std::array<T, 9> matmul(const std::array<T, 3>& lhs, const std::array<T, 3>& rhs) {
    return {
        lhs[0] * rhs[0], lhs[0] * rhs[1], lhs[0] * rhs[2],
        lhs[1] * rhs[0], lhs[1] * rhs[1], lhs[1] * rhs[2],
        lhs[2] * rhs[0], lhs[2] * rhs[1], lhs[2] * rhs[2]
    };
}

template<typename T>
std::array<T, 9> matmul(const std::array<T, 9>& lhs, const std::array<T, 3>& rhs) {
    return {
        lhs[0] * rhs[0], lhs[1] * rhs[1], lhs[2] * rhs[2],
        lhs[3] * rhs[0], lhs[4] * rhs[1], lhs[5] * rhs[2],
        lhs[6] * rhs[0], lhs[7] * rhs[1], lhs[8] * rhs[2]
    };
}

template<typename T>
std::array<T, 3> matmul(const std::array<T, 3>& vec, const std::array<T, 9>& mat) {
    return {
        vec[0] * mat[0] + vec[1] * mat[3] + vec[2] * mat[6], 
        vec[0] * mat[1] + vec[1] * mat[4] + vec[2] * mat[7],
        vec[0] * mat[2] + vec[1] * mat[5] + vec[2] * mat[8] 
    };
}

template<typename T>
std::array<T, 9> matmul(const std::array<T, 9>& a, const std::array<T, 9>& b) {
    return {
        // First row
        a[0]*b[0] + a[1]*b[3] + a[2]*b[6],
        a[0]*b[1] + a[1]*b[4] + a[2]*b[7],
        a[0]*b[2] + a[1]*b[5] + a[2]*b[8],
        
        // Second row
        a[3]*b[0] + a[4]*b[3] + a[5]*b[6],
        a[3]*b[1] + a[4]*b[4] + a[5]*b[7],
        a[3]*b[2] + a[4]*b[5] + a[5]*b[8],
        
        // Third row
        a[6]*b[0] + a[7]*b[3] + a[8]*b[6],
        a[6]*b[1] + a[7]*b[4] + a[8]*b[7],
        a[6]*b[2] + a[7]*b[5] + a[8]*b[8]
    };
}

template<typename T>
std::array<T, 3> vecmul(const std::array<T, 3>& lhs, const std::array<T, 3>& rhs) {
    return {
        lhs[1] * rhs[2] - lhs[2] * rhs[1], 
        rhs[0] * lhs[2] - lhs[0] * rhs[2],
        lhs[0] * rhs[1] - lhs[1] * rhs[0]
    };
}

template<typename T, size_t N>
std::array<T, N> operator*(T val, const std::array<T, N> arr) {
    return arr * val;
}

template<typename T>
constexpr std::array<T, 9> transpose(const std::array<T, 9>& matrix) {
    return {
        matrix[0], matrix[3], matrix[6], 
        matrix[1], matrix[4], matrix[7], 
        matrix[2], matrix[5], matrix[8] 
    };
}

template<typename T>
constexpr T div(const std::array<T, 9>& matrix) {
    return matrix[0] + matrix[4] + matrix[8];
}

}
