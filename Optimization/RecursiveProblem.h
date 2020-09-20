//
// Created by janos on 7/16/20.
//

#pragma once

#include "Functional.h"
#include "Tree.h"

namespace Phasefield::Solver {

struct RecursiveProblem {

    explicit RecursiveProblem(Tree& t);

    [[nodiscard]] size_t numParameters() const;

    [[nodiscard]] size_t numConstraints() const;

    void operator()(
            Containers::ArrayView<const double> parameters,
            double& cost,
            Containers::ArrayView<double> gradient,
            Containers::ArrayView<double> constraints,
            SparseMatrix* jacobian) const;

    void operator()(
            Containers::ArrayView<const double> parameters,
            double& cost,
            Containers::ArrayView<double> gradient) const;

    void determineSparsityStructure(SparseMatrix& jacobian) const;

    Tree& tree;
    size_t levelToOptimize = 0;

    Containers::Array<Functional> objectives;
    Containers::Array<Functional> constraints;
};

}

