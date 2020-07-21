//
// Created by janos on 7/16/20.
//

#pragma once

#include "modica_mortola.hpp"
#include "functional.h"

struct RecursiveFunctional {

    RecursiveFunctional(Containers::Array<const Mg::Vector3d> const& vs,
                        Containers::Array<const Mg::UnsignedInt> const& ts,
                        PhasefieldTree& t);

    Mg::UnsignedInt numParameters() const;
    Mg::UnsignedInt numConstraints() const;
    auto triangles() { return arrayCast<const Mg::Vector3ui>(indices); }

    void updateInternalDataStructures();

    void operator()(
            Containers::ArrayView<const double> data,
            double& cost,
            Containers::ArrayView<double> gradData,
            Containers::ArrayView<double> constr,
            SparseMatrix* jacobian) const;

    void getSparsityStructure(SparseMatrix& jacobian) const;

   Containers::Array<const Mg::Vector3d> const& vertices;
    Containers::Array<const Mg::UnsignedInt> const& indices;
    PhasefieldTree& tree;
    Containers::Array<double> weights;

    Containers::Array<Functional> objectives;
    Containers::Array<Functional> constraints;
};

