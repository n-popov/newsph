#pragma once

#include <string>
#include <fstream>
#include <nlohmann/json.hpp>

namespace config {

struct MaterialProperties {
    double density;
    double sound_speed_coefficient;
    double hugoniot_slope;
    double gruneisen_gamma;
    double shear_modulus;
    double yield_strength;
    double poisson_ratio;
    double bulk_modulus;
    double sound_speed;
    
    void calculate_derived_properties() {
        bulk_modulus = density * sound_speed_coefficient * sound_speed_coefficient;
        sound_speed = std::sqrt(bulk_modulus / density);
    }
};

struct SimulationParameters {
    double max_time;
    double dt;
    int output_frequency;
};

struct SPHParameters {
    double hdx;
    double dx;
    double h;
    double avisc_alpha;
    double avisc_beta;
    double avisc_eta;
    double xsph_eps;
    
    void calculate_derived_properties() {
        h = dx * hdx;
    }
};

struct ProjectileParameters {
    double radius;
    double velocity;
    std::string material;
};

struct PlateParameters {
    double length;
    double width;
    double thickness;
    std::string material;
};

struct GMSHParameters {
    bool is_enabled;
    bool verbose;
    bool write_mesh_file;
};

class SimulationConfig {
private:
    nlohmann::json json_data;
    
public:
    MaterialProperties aluminum_props;
    MaterialProperties steel_props;
    SimulationParameters simulation_params;
    SPHParameters sph_params;
    ProjectileParameters projectile_params;
    PlateParameters plate_params;
    GMSHParameters gmsh_params;
    
    void load_from_file(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Could not open config file: " + filename);
        }
        file >> json_data;
        file.close();
        
        parse_materials();
        parse_simulation_parameters();
        parse_sph_parameters();
        parse_projectile_parameters();
        parse_plate_parameters();
        parse_gmsh_parameters();
    }
    
    const MaterialProperties& get_material_properties(const std::string& material_name) const {
        if (material_name == "aluminum") {
            return aluminum_props;
        } else if (material_name == "steel") {
            return steel_props;
        } else {
            throw std::runtime_error("Unknown material: " + material_name);
        }
    }
    
private:
    void parse_materials() {
        auto& al = json_data["materials"]["aluminum"];
        aluminum_props.density = al["density"];
        aluminum_props.sound_speed_coefficient = al["sound_speed_coefficient"];
        aluminum_props.hugoniot_slope = al["hugoniot_slope"];
        aluminum_props.gruneisen_gamma = al["gruneisen_gamma"];
        aluminum_props.shear_modulus = al["shear_modulus"];
        aluminum_props.yield_strength = al["yield_strength"];
        aluminum_props.poisson_ratio = al["poisson_ratio"];
        aluminum_props.calculate_derived_properties();
        
        auto& st = json_data["materials"]["steel"];
        steel_props.density = st["density"];
        steel_props.sound_speed_coefficient = st["sound_speed_coefficient"];
        steel_props.hugoniot_slope = st["hugoniot_slope"];
        steel_props.gruneisen_gamma = st["gruneisen_gamma"];
        steel_props.shear_modulus = st["shear_modulus"];
        steel_props.yield_strength = st["yield_strength"];
        steel_props.poisson_ratio = st["poisson_ratio"];
        steel_props.calculate_derived_properties();
    }
    
    void parse_simulation_parameters() {
        auto& sim = json_data["simulation"];
        simulation_params.max_time = sim["max_time"];
        simulation_params.dt = sim["dt"];
        simulation_params.output_frequency = sim["output_frequency"];
    }
    
    void parse_sph_parameters() {
        auto& sph = json_data["sph_parameters"];
        sph_params.hdx = sph["hdx"];
        sph_params.dx = sph["dx"];
        sph_params.avisc_alpha = sph["avisc_alpha"];
        sph_params.avisc_beta = sph["avisc_beta"];
        sph_params.avisc_eta = sph["avisc_eta"];
        sph_params.xsph_eps = sph["xsph_eps"];
        sph_params.calculate_derived_properties();
    }
    
    void parse_projectile_parameters() {
        auto& proj = json_data["projectile"];
        projectile_params.radius = proj["radius"];
        projectile_params.velocity = proj["velocity"];
        projectile_params.material = proj["material"];
    }
    
    void parse_plate_parameters() {
        auto& plate = json_data["plate"];
        plate_params.length = plate["length"];
        plate_params.width = plate["width"];
        plate_params.thickness = plate["thickness"];
        plate_params.material = plate["material"];
    }

    void parse_gmsh_parameters() {
        auto& gmsh = json_data["gmsh"];

        gmsh_params.is_enabled = gmsh["is_enabled"];
        gmsh_params.verbose = gmsh["verbose"];
        gmsh_params.write_mesh_file = gmsh["write_mesh_file"];
    }
};

} // namespace config
