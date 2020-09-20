//
// Created by janos on 27.11.19.
//

#pragma once

#include "Types.h"
#include "Surface.h"
#include "Functional.h"

namespace Phasefield {

struct DirichletEnergy {

    explicit DirichletEnergy(Mesh& mesh);

    [[nodiscard]] size_t numParameters() const;

    [[nodiscard]] static FunctionalType::Value type() { return FunctionalType::DirichletEnergy; }

    template<class Scalar>
    void operator()(ArrayView<const Scalar> const& parameters,
                    ArrayView<const Scalar> const& weights,
                    Scalar& out,
                    ArrayView<Scalar> const& gradP,
                    ArrayView<Scalar> const& gradW);

    Mesh& mesh;
};

struct AreaRegularizer {

    explicit AreaRegularizer(Mesh& mesh);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> const& parameters,
                    ArrayView<const Scalar> const& weights,
                    Scalar& out,
                    ArrayView<Scalar> const& gradP,
                    ArrayView<Scalar> const& gradW);

    [[nodiscard]] size_t numParameters() const;

    [[nodiscard]] static FunctionalType::Value type() { return FunctionalType::AreaRegularizer; }

    double areaRatio = 0.5;
    Mesh& mesh;
};

struct DoubleWellPotential {

    explicit DoubleWellPotential(Mesh& mesh);

    template<class Scalar>
    void operator()(ArrayView<const Scalar> const& parameters,
                    ArrayView<const Scalar> const& weights,
                    Scalar& out,
                    ArrayView<Scalar> const& gradP,
                    ArrayView<Scalar> const& gradW);

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

}


