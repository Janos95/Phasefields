//
// Created by janos on 8/25/20.
//

#pragma once

#include "Surface.h"
#include "Types.h"
#include "Functional.h"

namespace Phasefield {

struct WeakYamabe {

    explicit WeakYamabe(Mesh& m);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> const& parameters,
                    ArrayView<const Scalar> const& weights,
                    double& out,
                    ArrayView<Scalar> const& gradP,
                    ArrayView<Scalar> const& gradW);

    [[nodiscard]] size_t numParameters() const;

    Mesh* mesh;
};

DECLARE_FUNCTIONAL_CONSTRUCTOR(WeakYamabe)
DECLARE_FUNCTIONAL_OPERATOR(WeakYamabe, double)

}

