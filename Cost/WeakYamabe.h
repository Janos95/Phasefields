//
// Created by janos on 8/25/20.
//

#pragma once

#include "Surface.h"
#include "Types.h"
#include "Functional.h"

class adouble;

namespace Phasefield {

struct WeakYamabe {

    explicit WeakYamabe(Mesh& m);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> parameters,
                    ArrayView<const Scalar> weights,
                    Scalar& out,
                    ArrayView<Scalar> gradP,
                    ArrayView<Scalar> gradW);

    [[nodiscard]] size_t numParameters() const;

    Mesh* mesh;
};

DECLARE_FUNCTIONAL_CONSTRUCTOR(WeakYamabe)
DECLARE_FUNCTIONAL_OPERATOR(WeakYamabe, double)
DECLARE_FUNCTIONAL_OPERATOR(WeakYamabe, adouble)

}

