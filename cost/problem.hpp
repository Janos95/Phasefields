//
// Created by janos on 28.03.20.
//

#pragma once

#include "unique_function.h"
#include "functional.hpp"

#include <Magnum/Trade/Trade.h>
#include <Corrade/Containers/Containers.h>
#include <mutex>

namespace Cr = Corrade;

namespace solver {

struct Problem {

    //hack these in :/
    VisualizationFlags flags = VisualizationFlag::Phasefield;
    VisualizationFlags* update = nullptr;
    std::mutex* mutex = nullptr;
    Mg::Trade::MeshData* meshData;

    [[nodiscard]] int numParameters() const;
    [[nodiscard]] int numConstraints() const;


    Mg::Double h = 1e-5;
    void numericalGradient(double* parameters, double *jacF, double* jacC);

    Cr::Containers::Array<Cr::Containers::Pointer<FunctionalD>> functionals;
    Cr::Containers::Array<Cr::Containers::Pointer<FunctionalD>> constraints;

    bool evaluate(
            const double *parameters,
            double *costFunctional,
            double *jacobiansFunctionals,
            double *costConstraint,
            double *jacobiansConstraints) const;

    bool evaluateObjective(const double *parameters, double *cost, double *jacobians) const;

    bool evaluateConstraints(const double *parameters, double *cost, double *jacobians) const;

};

}