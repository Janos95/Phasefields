//
// Created by janos on 28.03.20.
//

#pragma once


#include "visualization_proxy.hpp"
#include "functional.h"
#include "optimization.h"

#include <Corrade/Containers/Array.h>


namespace solver {

struct Problem {

    [[nodiscard]] int numParameters() const;
    [[nodiscard]] int numConstraints() const;

    Containers::Array<Functional> objectives;
    Containers::Array<Functional> constraints;

    VisualizationProxy proxy;

    bool evaluateObjective(const double *parameters, double *cost, SparseMatrix* gradient, SparseMatrix* hessian) const;

    bool evaluateConstraints(const double *parameters, double *cost, SparseMatrix* gradient, SparseMatrix* hessian, double* scaling) const;

};

}