//
// Created by janos on 27.11.19.
//

#pragma once

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/StridedArrayView.h>
#include "c1_functions.hpp"
#include "Enums.h"

namespace Phasefield {

struct DirichletEnergy {
    DirichletEnergy(
            Containers::Array<Mg::Vector3d> const& vertices,
            Containers::Array<Mg::UnsignedInt> const& indices);

    uint32_t numParameters() const;

    template<class Scalar>
    void operator()(Containers::ArrayView<const Scalar> const& parameters,
                    Containers::ArrayView<const Scalar> const& weights,
                    Scalar& out,
                    Containers::ArrayView<Scalar> const& gradP,
                    Containers::ArrayView<Scalar> const& gradW);

    void updateInternalDataStructures();

    Containers::Array<const Mg::UnsignedInt> const& indices;
    Containers::Array<const Mg::Vector3d> const& vertices;

    Containers::Array<Mg::Double> areas;

    Containers::Array<Mg::Vector3d> edgeNormalsData;
    Containers::StridedArrayView2D<Mg::Vector3d> edgeNormals;
};

struct AreaRegularizer {

    AreaRegularizer(Containers::Array<Magnum::Vector3d> const&, Containers::Array<Mg::UnsignedInt> const&);

    void updateInternalDataStructures();

    template<class Scalar>
    void operator()(Containers::ArrayView<const Scalar> const& parameters,
                    Containers::ArrayView<const Scalar> const& weights,
                    Scalar& out,
                    Containers::ArrayView<Scalar> const& gradP,
                    Containers::ArrayView<Scalar> const& gradW);

    mutable Mg::Double currentArea;
    Mg::Double area, areaRatio = 0.5;

    Containers::Array<Mg::UnsignedInt> const& indices;
    Containers::Array<Mg::Vector3d> const& vertices;
    Containers::Array<Mg::Double> integralOperator;
};

struct DoubleWellPotential {

    DoubleWellPotential(Containers::Array<Magnum::Vector3d> const&,
                        Containers::Array<Mg::UnsignedInt> const&);

    void updateInternalDataStructures();

    template<class Scalar>
    void operator()(Containers::ArrayView<const Scalar> const& parameters,
                    Containers::ArrayView<const Scalar> const& weights,
                    Scalar& out,
                    Containers::ArrayView<Scalar> const& gradP,
                    Containers::ArrayView<Scalar> const& gradW);

    Containers::Array<Mg::UnsignedInt> const& indices;
    Containers::Array<Mg::Vector3d> const& vertices;
    Containers::Array<Mg::Double> integralOperator;
};


extern template void DirichletEnergy::operator()(Containers::ArrayView<const double> const& parameters,
                                                 Containers::ArrayView<const double> const& weights,
                                                 double& out,
                                                 Containers::ArrayView<double> const& gradP,
                                                 Containers::ArrayView<double> const& gradW);

extern template void AreaRegularizer::operator()(Containers::ArrayView<const double> const& parameters,
                                                 Containers::ArrayView<const double> const& weights,
                                                 double& out,
                                                 Containers::ArrayView<double> const& gradP,
                                                 Containers::ArrayView<double> const& gradW);

extern template void DoubleWellPotential::operator()(Containers::ArrayView<const double> const& parameters,
                                                     Containers::ArrayView<const double> const& weights,
                                                     double& out,
                                                     Containers::ArrayView<double> const& gradP,
                                                     Containers::ArrayView<double> const& gradW);

}
