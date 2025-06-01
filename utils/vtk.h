#pragma once

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>
#include <filesystem>

#include "../method/particle.h"

using namespace mysph;

template <typename T>
void write_particles_vtk(const std::string& filename, const std::vector<Particle<T>>& particles) {
    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble();
    for (const auto& p : particles) {
        points->InsertNextPoint(p.r[0], p.r[1], p.r[2]);
    }

    auto vertices = vtkSmartPointer<vtkCellArray>::New();
    for (vtkIdType i = 0; i < particles.size(); ++i) {
        vertices->InsertNextCell(1, &i);
    }

    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetVerts(vertices);

    auto add_vector_field = [&](const char* name, const vec3<T> Particle<T>::*field) {
        auto array = vtkSmartPointer<vtkFloatArray>::New();
        array->SetName(name);
        array->SetNumberOfComponents(3);
        array->SetNumberOfTuples(particles.size());

        for (size_t i = 0; i < particles.size(); ++i) {
            const auto& val = particles[i].*field;
            array->SetTuple3(i, val[0], val[1], val[2]);
        }
        polyData->GetPointData()->AddArray(array);
    };

    auto add_scalar_field = [&](const char* name, const T Particle<T>::*field) {
        auto array = vtkSmartPointer<vtkFloatArray>::New();
        array->SetName(name);
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(particles.size());

        for (size_t i = 0; i < particles.size(); ++i) {
            array->SetTuple1(i, particles[i].*field);
        }
        polyData->GetPointData()->AddArray(array);
    };

    auto add_int_field = [&](const char* name, const int Particle<T>::*field) {
        auto array = vtkSmartPointer<vtkFloatArray>::New();
        array->SetName(name);
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(particles.size());

        for (size_t i = 0; i < particles.size(); ++i) {
            array->SetTuple1(i, static_cast<float>(particles[i].*field));
        }
        polyData->GetPointData()->AddArray(array);
    };

    add_vector_field("velocity", &Particle<T>::v);
    add_vector_field("vstar", &Particle<T>::vstar);
    add_vector_field("Fv", &Particle<T>::Fv);
    add_vector_field("F", &Particle<T>::F);

    add_scalar_field("density", &Particle<T>::rho);
    add_scalar_field("reference_density", &Particle<T>::rho0);
    add_scalar_field("pressure", &Particle<T>::p);
    add_scalar_field("mass", &Particle<T>::m);
    add_scalar_field("soundspeed", &Particle<T>::cs);
    add_scalar_field("bulk_modulus", &Particle<T>::k);
    add_scalar_field("shear_modulus", &Particle<T>::G);
    add_scalar_field("yield_stress", &Particle<T>::Yo);
    add_scalar_field("J2", &Particle<T>::J2);
    add_scalar_field("plastic_strain", &Particle<T>::plastic_strain);
    add_int_field("material", &Particle<T>::material);

    // Add stress tensor components
    add_scalar_field("stress_00", &Particle<T>::s00);
    add_scalar_field("stress_01", &Particle<T>::s01);
    add_scalar_field("stress_02", &Particle<T>::s02);
    add_scalar_field("stress_11", &Particle<T>::s11);
    add_scalar_field("stress_12", &Particle<T>::s12);
    add_scalar_field("stress_22", &Particle<T>::s22);

    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(polyData);
    writer->Write();
}
