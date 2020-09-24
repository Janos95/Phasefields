//
// Created by janos on 27.11.19.
//

#pragma once

#include "Types.h"
#include "Surface.h"
#include "Functional.h"

class adouble;

namespace Phasefield {

struct DirichletEnergy {

    explicit DirichletEnergy(Mesh& mesh);

    [[nodiscard]] size_t numParameters() const;

    [[nodiscard]] static FunctionalType::Value type() { return FunctionalType::DirichletEnergy; }

    template<class Scalar>
    void operator()(ArrayView<const Scalar> parameters,
                    ArrayView<const Scalar> weights,
                    Scalar& out,
                    ArrayView<Scalar> gradP,
                    ArrayView<Scalar> gradW);

    Mesh& mesh;
};

struct AreaRegularizer {

    explicit AreaRegularizer(Mesh& mesh);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> parameters,
                    ArrayView<const Scalar> weights,
                    Scalar& out,
                    ArrayView<Scalar> gradP,
                    ArrayView<Scalar> gradW);

    [[nodiscard]] size_t numParameters() const;

    [[nodiscard]] static FunctionalType::Value type() { return FunctionalType::AreaRegularizer; }

    double targetArea = 1;
    Mesh& mesh;
};

struct DoubleWellPotential {

    explicit DoubleWellPotential(Mesh& mesh);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> parameters,
                    ArrayView<const Scalar> weights,
                    Scalar& out,
                    ArrayView<Scalar> gradP,
                    ArrayView<Scalar> gradW);

    [[nodiscard]] size_t numParameters() const;

    [[nodiscard]] static FunctionalType::Value type() { return FunctionalType::DoubleWellPotential; }

    Mesh& mesh;
};

DECLARE_FUNCTIONAL_CONSTRUCTOR(DirichletEnergy)
DECLARE_FUNCTIONAL_CONSTRUCTOR(AreaRegularizer)
DECLARE_FUNCTIONAL_CONSTRUCTOR(DoubleWellPotential)

DECLARE_FUNCTIONAL_OPERATOR(DirichletEnergy, double)
DECLARE_FUNCTIONAL_OPERATOR(AreaRegularizer, double)
DECLARE_FUNCTIONAL_OPERATOR(DoubleWellPotential, double)

DECLARE_FUNCTIONAL_OPERATOR(DirichletEnergy, adouble)
DECLARE_FUNCTIONAL_OPERATOR(AreaRegularizer, adouble)
DECLARE_FUNCTIONAL_OPERATOR(DoubleWellPotential, adouble)


}


