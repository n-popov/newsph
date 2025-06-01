cmake_minimum_required(VERSION 3.10)
project(SPH_Impact)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/opt/homebrew/lib/libgmsh.dylib"
    INTERFACE_INCLUDE_DIRECTORIES "/usr/local/include"
)
set(GMSH_LIBRARIES gmsh)

# Find VTK
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include(FetchContent)
FetchContent_Declare(json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG v3.11.2
)

# Create your executable
add_executable(sph_impact main.cpp)
target_link_libraries(sph_impact ${VTK_LIBRARIES} ${GMSH_LIBRARIES})

# Copy config file to build directory
configure_file(${CMAKE_SOURCE_DIR}/config.json ${CMAKE_BINARY_DIR}/config.json COPYONLY)
