//
// Created by janos on 7/16/20.
//

#pragma once

#include "Functional.h"
#include "PhasefieldTree.h"

namespace Solver {

struct RecursiveProblem {

    explicit RecursiveProblem(PhasefieldTree& t);

    Mg::UnsignedInt numParameters() const;

    Mg::UnsignedInt numConstraints() const;

    void updateInternalDataStructures();

    void operator()(
            Containers::ArrayView<const double> parameters,
            double& cost,
            Containers::ArrayView<double> gradient,
            Containers::ArrayView<double> constraints,
            SparseMatrix* jacobian) const;

    void determineSparsityStructure(SparseMatrix& jacobian) const;

    PhasefieldTree& tree;

    Containers::Array<Functional> objectives;
    Containers::Array<Functional> constraints;
};

}

