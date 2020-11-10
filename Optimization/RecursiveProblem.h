//
// Created by janos on 7/16/20.
//

#pragma once

#include "Functional.h"
#include "Tree.h"

#include <Corrade/Containers/StridedArrayView.h>

namespace Phasefield::Solver {

struct RecursiveProblem {

    explicit RecursiveProblem(Tree& t);

    [[nodiscard]] size_t numParameters() const;

    //[[nodiscard]] size_t numConstraints() const;

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

    //void determineSparsityStructure(SparseMatrix& jacobian) const;

    Tree& tree;
    Node nodeToOptimize{Invalid, nullptr};

    struct Objective {
        Functional f{nullptr};
        Array<Vector2> history;
        bool show = true;
    };

    bool evaluateObjective = true;

    Array<Objective> objectives;
    Array<Functional> constraints;
};

}

