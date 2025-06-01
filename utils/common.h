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

template<typename T, size_t N>
std::array<T, N> operator*(T val, const std::array<T, N> arr) {
    return arr * val;
}

template<typename T>
constexpr std::array<T, 9> transpose(const std::array<T, 9>& matrix) {
    return {
        matrix[0], matrix[3], matrix[6],  // First column becomes first row
        matrix[1], matrix[4], matrix[7],  // Second column becomes second row
        matrix[2], matrix[5], matrix[8]   // Third column becomes third row
    };
}

template<typename T>
constexpr T div(const std::array<T, 9>& matrix) {
    return matrix[0] + matrix[4] + matrix[8];
}


}
