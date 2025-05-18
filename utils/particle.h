// utils/particle.h
#pragma once

#include <array>
#include <numeric>
#include <vector>

namespace mysph {

template<typename T>
using vec3 = std::array<T, 3>;

template<typename T>
struct Particle {
    Particle() = default;
    Particle(vec3<T> r, vec3<T> v, T m, T rho, int material, T k, T nu)
        : r(r), v(v), m(m), rho(rho), material(material), k(k), nu(nu) {}
    Particle(vec3<T> r) : r(r) {}

    // Positions and velocities
    vec3<T> r = {T(0), T(0), T(0)};
    vec3<T> v = {T(0), T(0), T(0)};
    vec3<T> vstar = {T(0), T(0), T(0)};
    
    // Forces
    vec3<T> Fv = {T(0), T(0), T(0)};
    vec3<T> F = {T(0), T(0), T(0)};
    
    // Basic properties
    T m = T(0);          // mass
    T rho = T(0);        // density
    T rho0 = T(0);       // reference density
    T p = T(0);          // pressure
    T e = T(0);          // internal energy
    T cs = T(0);         // speed of sound
    
    // Material properties
    int material = 0;    // material type (0=plate, 1=projectile)
    T k = T(0);          // bulk modulus
    T nu = T(0);         // viscosity
    T G = T(0);          // shear modulus
    T Yo = T(0);         // yield stress
    
    // Stress tensor components (3D)
    T s00 = T(0), s01 = T(0), s02 = T(0);
    T s10 = T(0), s11 = T(0), s12 = T(0);
    T s20 = T(0), s21 = T(0), s22 = T(0);
    
    // Velocity gradient tensor
    T v00 = T(0), v01 = T(0), v02 = T(0);
    T v10 = T(0), v11 = T(0), v12 = T(0);
    T v20 = T(0), v21 = T(0), v22 = T(0);
    
    // Artificial stress components
    T as00 = T(0), as01 = T(0), as02 = T(0);
    T as11 = T(0), as12 = T(0), as22 = T(0);
    
    // For plasticity calculations
    T plastic_strain = T(0);
    T J2 = T(0);         // Second invariant of stress deviator
};

}
