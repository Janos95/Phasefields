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
    void operator()(Scalar const*, Scalar*, SparseMatrix*) const;

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

    uint32_t numParameters() const  {
        return vertices.size();
    }

    void updateInternalDataStructures()  {
        integralOperator = computeIntegralOperator(triangles(), vertices);
    }

    template<class Scalar>
    void operator()(Scalar const* params,
                  Scalar* cost,
                  SparseMatrix* jacobians,
                  SparseMatrix* hessian) const {

        *cost = 0.;
        for (int i = 0; i < vertices.size(); ++i){
            *cost += f.eval(params[i]) * integralOperator[i];
            if(jacobians)
                jacobians->values[i] = f.grad(params[i]) * integralOperator[i];
            if(hessian)
                hessian->values[i] = f.hess(params[i]) * integralOperator[i];
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
    void operator()(Scalar const* params, Scalar* cost, SparseMatrix* jacobians, SparseMatrix* hessian) const;

    IntegralFunctional<SmootherStep> integral;
    mutable Mg::Double currentArea;
    Mg::Double area, areaRatio = 0.5;
};

struct AreaConstraints {

    AreaConstraints(Containers::Array<const Mg::Vector3d> const& vs,
                    Containers::Array<const Mg::UnsignedInt> const& ts,
                    PhasefieldTree& t);

    Mg::UnsignedInt numParameters() const;
    Mg::UnsignedInt numResiduals() const;
    auto triangles() { return arrayCast<const Mg::Vector3ui>(indices); }

    void updateInternalDataStructures();

    template<class Scalar>
    void operator()(Scalar const* params, Scalar* cost) const;

    PhasefieldTree& tree;
    Containers::Array<double> weights;
    Containers::Array<const Mg::Vector3d> const& vertices;
    Containers::Array<const Mg::UnsignedInt> const& indices;
};

using DoubleWellPotential = IntegralFunctional<DoubleWell>;

class adouble;

extern template void DirichletEnergy::operator()(adouble const*, adouble*, SparseMatrix*) const;
extern template void DirichletEnergy::operator()(double const*, double*, SparseMatrix*) const;

extern template void AreaRegularizer::operator()(adouble const*, adouble*, SparseMatrix*, SparseMatrix*) const;
extern template void AreaRegularizer::operator()(double const*, double*, SparseMatrix*, SparseMatrix*) const;

extern template void AreaConstraints::operator()(adouble const*, adouble*) const;
extern template void AreaConstraints::operator()(double const*, double*) const;

extern template void DoubleWellPotential::operator()(adouble const*, adouble*, SparseMatrix*, SparseMatrix*) const;
extern template void DoubleWellPotential::operator()(double const*, double*, SparseMatrix*, SparseMatrix*) const;


