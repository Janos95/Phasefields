//
// Created by janos on 2/23/20.
//
#include "modica_mortola.hpp"

#include <Corrade/Utility/Assert.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/FunctionsBatch.h>

using namespace Magnum;
using namespace Corrade;

template<class Scalar>
DirichletEnergy<Scalar>::DirichletEnergy(
        Containers::Array<const Vector3d> const& vs,
        Containers::Array<const UnsignedInt> const& is):
    indices(is),
    vertices(vs),
    areas(computeAreas(triangles(), vertices)),
    stiffnessMatrix(computeStiffnessMatrix(triangles(), vs)),
    diagonal(stiffnessMatrix.diagonal())
{
    CORRADE_ASSERT(!Math::isNan(areas), "Dirichlet Energy : areas contains NaN",);
}


template<class Scalar>
uint32_t DirichletEnergy<Scalar>::numParameters() const {
    return vertices.size();
}

template<class Scalar>
bool DirichletEnergy<Scalar>::evaluate(double const* parameters,
              double* residual,
              double* jacobians) const
{
    Eigen::Map<const Eigen::VectorXd> phasefield(parameters, numParameters());
    auto halfGrad = stiffnessMatrix * phasefield;
    if(jacobians){
        Eigen::Map<Eigen::VectorXd> jac(jacobians, numParameters());
        halfGrad.eval();
        jac = 2 * halfGrad;
    }
    *residual = phasefield.transpose() * halfGrad;
    return true;
}


