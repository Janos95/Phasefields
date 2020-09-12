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


extern template void DirichletEnergy::operator()(ArrayView<const double> const& parameters,
                                                 ArrayView<const double> const& weights,
                                                 double& out,
                                                 ArrayView<double> const& gradP,
                                                 ArrayView<double> const& gradW);

extern template void AreaRegularizer::operator()(ArrayView<const double> const& parameters,
                                                 ArrayView<const double> const& weights,
                                                 double& out,
                                                 ArrayView<double> const& gradP,
                                                 ArrayView<double> const& gradW);

extern template void DoubleWellPotential::operator()(ArrayView<const double> const& parameters,
                                                     ArrayView<const double> const& weights,
                                                     double& out,
                                                     ArrayView<double> const& gradP,
                                                     ArrayView<double> const& gradW);


extern template Functional::Functional(DirichletEnergy);
extern template Functional::Functional(AreaRegularizer);
extern template Functional::Functional(DoubleWellPotential);

}


