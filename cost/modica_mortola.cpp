//
// Created by janos on 2/23/20.
//
#include "modica_mortola.hpp"
#include "interpolation.hpp"

#include <Corrade/Utility/Assert.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/FunctionsBatch.h>

using namespace Magnum;
using namespace Corrade;

DirichletEnergy::DirichletEnergy(
        Containers::ArrayView<const Vector3d> const& vertices_,
        Containers::ArrayView<const Vector3ui> const& triangles_):
    Functional(Functional::MetaData::AllocateFromLoss(TrivialLoss{}), FunctionalType::DirichletEnergy),
    triangles(triangles_),
    vertices(vertices_),
    areas(computeAreas(triangles_, vertices_)),
    stiffnessMatrix(computeStiffnessMatrix(triangles_, vertices_)),
    diagonal(stiffnessMatrix.diagonal())
{
    CORRADE_ASSERT(!Math::isNan(areas), "Dirichlet Energy : areas contains NaN",);
}


int DirichletEnergy::numParameters() const {
    return vertices.size();
}

bool DirichletEnergy::evaluate(double const* parameters,
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


