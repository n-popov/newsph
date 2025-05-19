#pragma once

#include <array>
#include <cmath>
#include <numeric>

namespace mysph {

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

}
