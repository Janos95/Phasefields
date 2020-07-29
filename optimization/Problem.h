//
// Created by janos on 28.03.20.
//

#pragma once


#include "visualization_proxy.hpp"
#include "Functional.h"
#include "Optimization.h"
#include "SparseMatrix.h"

#include <Corrade/Containers/Array.h>

namespace Solver {

struct Problem {

    explicit Problem();

    ~Problem();

    [[nodiscard]] std::size_t numParameters() const;

    [[nodiscard]] std::size_t numConstraints() const;

    Containers::Array<Functional> objectives, constraints;
    SparseMatrix hessian, jacobian; /* hessian of the lagrangian, jacobian of the constraints */
    VisualizationProxy* proxy = nullptr; /* pipe to visualize stuff */
    int tagL, tagJ, tagG;

    void updateInternalDataStructures(const Containers::Array<double>&);

    void evaluate(
            double const* parameters,
            double* residual,
            double* constr,
            double* g,
            SparseMatrix* j,
            SparseMatrix* h,
            double objectiveScale,
            double const* lambdas) const;
};

}