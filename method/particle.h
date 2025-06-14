#pragma once

#include <array>
#include <numeric>
#include <vector>

namespace mysph {

template<typename T>
using vec3 = std::array<T, 3>;

template<typename T>
using vec9 = std::array<T, 9>;

template<typename T>
struct Particle {
    Particle() = default;
    Particle(vec3<T> r, vec3<T> v, T m, T rho, int material, T k, T nu)
        : r(r), v(v), m(m), rho(rho), material(material), k(k), nu(nu) {}
    Particle(vec3<T> r) : r(r) {}

    std::vector<Particle<T>*>& get_neighbors() {
        if (!neighbors.empty()) return neighbors;

        if (external_neighbors.has_value()) return *(external_neighbors.value());

        throw std::runtime_error("either neighbors or external_neighbors should be specified");
    }

    // Positions and velocities
    vec3<T> r = {T(0), T(0), T(0)};
    vec3<T> v = {T(0), T(0), T(0)};
    vec3<T> vstar = {T(0), T(0), T(0)};
    
    // Forces
    vec3<T> Fv = {T(0), T(0), T(0)};
    
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
    vec9<T> stress = {};
    
    // Velocity gradient tensor
    vec9<T> v_grad = {};
    
    // Artificial stress components
    vec9<T> a_stress = {};
    
    // For plasticity calculations
    T plastic_strain = T(0);
    T J2 = T(0);         // Second invariant of stress deviator

    vec3<T> a = {}; // acceleration
    T ae = T(0); // energy change rate
    vec9<T> acc_stress = {}; // stress change rate
    T arho = T(0); // density change rate

    // Fake particles
    bool is_fake = false;

    double cf = 0.;

    std::vector<Particle<T>*> neighbors;
    std::optional<std::vector<Particle<T>*>*> external_neighbors;
};

}
