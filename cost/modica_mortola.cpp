//
// Created by janos on 2/23/20.
//
#include "modica_mortola.hpp"
#include "sparse_matrix.h"
#include "fem.hpp"

#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/Array.h>

#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <adolc/adouble.h>
#include <eigen3/unsupported/Eigen/AdolcForward>
#include <numeric>
#include <algorithm>

using namespace Magnum;
using namespace Corrade;


DirichletEnergy::DirichletEnergy(
        Containers::Array<const Vector3d> const& vs,
        Containers::Array<const UnsignedInt> const& is):
    indices(is),
    vertices(vs),
    areas(computeAreas(triangles(), vs)),
    stiffnessMatrix(computeStiffnessMatrix(triangles(), vs)),
    diagonal(stiffnessMatrix.diagonal())
{
    CORRADE_ASSERT(!Math::isNan(areas), "Dirichlet Energy : areas contains NaN",);
}


uint32_t DirichletEnergy::numParameters() const {
    return vertices.size();
}

adouble test;

template<class Scalar>
void DirichletEnergy::operator()(Scalar const* parameters, Scalar* residual, Scalar* gradient) const
{
    auto n = numParameters();
    Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> phasefield(parameters, n);
    auto halfGrad = stiffnessMatrix.cast<Scalar>() * phasefield;
    if(gradient){
        Eigen::Map<Eigen::VectorXd> jac(gradient, numParameters());
        halfGrad.eval();
        for(std::size_t i = 0; i < n; ++i){
            if constexpr(std::is_same_v<Scalar, adouble>)
                jac[i] = 2 * halfGrad[i].getValue();
            else
                jac[i] = 2 * halfGrad[i];
        }
    }
    *residual = phasefield.transpose() * halfGrad;
}

AreaRegularizer::AreaRegularizer(
        Containers::Array<const Magnum::Vector3d> const& vs,
        Containers::Array<const Mg::UnsignedInt> const& is) :
    integral(vs, is)
{
    auto areas = computeAreas(integral.triangles(), integral.vertices);
    area = std::accumulate(areas.begin(), areas.end(), 0.);
}

void AreaRegularizer::updateInternalDataStructures()  {
    integral.updateInternalDataStructures();
    auto areas = computeAreas(integral.triangles(), integral.vertices);
    area = std::accumulate(areas.begin(), areas.end(), 0.);
}

template<class Scalar>
void AreaRegularizer::operator()(Scalar const* params, Scalar* cost, SparseMatrix* jacobians, SparseMatrix* hessian) const {
    Scalar c = currentArea;
    integral(params, &c, jacobians, hessian);
    *cost = c - areaRatio * area; /*@todo racy */
}



template void DirichletEnergy::operator()(adouble const*, adouble*, SparseMatrix*) const;
template void DirichletEnergy::operator()(double const*, double*, SparseMatrix*) const;

template void AreaRegularizer::operator()(adouble const*, adouble*, SparseMatrix*, SparseMatrix*) const;
template void AreaRegularizer::operator()(double const*, double*, SparseMatrix*, SparseMatrix*) const;

template void DoubleWellPotential::operator()(adouble const*, adouble*, SparseMatrix*, SparseMatrix*) const;
template void DoubleWellPotential::operator()(double const*, double*, SparseMatrix*, SparseMatrix*) const;
