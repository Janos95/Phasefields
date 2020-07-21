//
// Created by janos on 27.11.19.
//

#pragma once

#include "c1_functions.hpp"
#include "fem.hpp"
#include "y_combinator.hpp"
#include "phasefield_tree.hpp"
#include "types.hpp"
#include "sparse_matrix.h"

#include <Eigen/SparseCore>

struct DirichletEnergy
{
    DirichletEnergy(
            Containers::Array<const Mg::Vector3d> const& vertices,
            Containers::Array<const Mg::UnsignedInt> const& indices);

    uint32_t numParameters() const;

    template<class Scalar>
    void operator()(Scalar const* params,
                    Scalar const* weights,
                    Scalar* cost) const;

        auto triangles() { return arrayCast<Mg::Vector3ui>(indices); }

    Containers::Array<const Mg::UnsignedInt> const& indices;
    Containers::Array<const Mg::Vector3d> const& vertices;
    Containers::Array<Mg::Double> areas;

    Eigen::SparseMatrix<double> stiffnessMatrix;
    Eigen::Matrix<double, Eigen::Dynamic, 1> diagonal;
};

template<class F>
struct IntegralFunctional
{
    IntegralFunctional(
            Containers::Array<const Mg::Vector3d> const& vs,
            Containers::Array<const Mg::UnsignedInt> const& is):
            indices(is),
            vertices(vs),
            integralOperator(computeIntegralOperator(triangles(), vertices))
    {
    }

    uint32_t numParameters() const {
        return vertices.size();
    }

    void updateInternalDataStructures()  {
        integralOperator = computeIntegralOperator(triangles(), vertices);
    }

    template<class Scalar>
    void operator()(Scalar const* params,
                    Scalar const* weights,
                    Scalar* cost) const {
        *cost = 0.;
        for (int i = 0; i < vertices.size(); ++i){
            *cost += f.eval(params[i]) * weights[i] * integralOperator[i];
        }
    }

    template<class Scalar>
    void operator()(Containers::ArrayView<const Scalar> parameters,
            Containers::ArrayView<const Scalar> weights,
            double& cost,
            Containers::ArrayView<Scalar> gradP,
            Containers::ArrayView<Scalar> gradW)
    {
        for (int i = 0; i < vertices.size(); ++i){
            cost += f.eval(parameters[i]) * weights[i] * integralOperator[i];
            if(gradP){
                gradP[i] += f.grad(parameters[i]) * weights[i] * integralOperator[i];
            }
            if(gradW){
                gradW[i] += f.eval(parameters[i]) * integralOperator[i];
            }
        }
    }

    auto triangles() { return arrayCast<Mg::Vector3ui>(indices); }

    F f;
    Containers::Array<const Mg::UnsignedInt> const& indices;
    Containers::Array<const Mg::Vector3d> const& vertices;
    Containers::Array<Mg::Double> integralOperator;
};

struct AreaRegularizer{

    AreaRegularizer(Containers::Array<const Magnum::Vector3d> const&, Containers::Array<const Mg::UnsignedInt> const&);

    void updateInternalDataStructures();

    template<class Scalar>
    void operator()(Containers::ArrayView<const Scalar> const& parameters,
            Containers::ArrayView<const Scalar> const& weights,
            Scalar& out,
            Containers::ArrayView<Scalar> const& gradP,
            Containers::ArrayView<Scalar> const& gradW);

    IntegralFunctional<SmootherStep> integral;

    mutable Mg::Double currentArea;
    Mg::Double area, areaRatio = 0.5;
};

using DoubleWellPotential = IntegralFunctional<DoubleWell>;

class adouble;

extern template void DirichletEnergy::operator()(adouble const*, adouble const*, adouble*) const;

extern template void AreaRegularizer::operator()(adouble const*,adouble const*, adouble*) const;

extern template void DoubleWellPotential::operator()(adouble const*, adouble const*, adouble*) const;


