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

    bool evaluate(const double *parameters,
                  double *cost,
                  double *jacobians) const;

    [[nodiscard]] int numParameters() const;

    Cr::Containers::Array<Cr::Containers::Pointer<Functional>> functionals;
    Cr::Containers::Array<Cr::Containers::Pointer<Functional>> constraints;
};

}