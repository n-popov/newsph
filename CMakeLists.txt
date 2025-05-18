cmake_minimum_required(VERSION 3.10)
project(SPH_Impact)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(sph_impact main.cpp)
target_link_libraries(sph_impact ${VTK_LIBRARIES})
