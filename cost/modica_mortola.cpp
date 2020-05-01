//
// Created by janos on 2/23/20.
//
#include "modica_mortola.hpp"
#include "interpolation.hpp"

#include <Corrade/Utility/Assert.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>


using namespace Magnum;
using namespace Corrade;

DirichletEnergy::DirichletEnergy(
        Containers::ArrayView<const Vector3d> const& vertices_,
        Containers::ArrayView<const Vector3ui> const& triangles_):
    triangles(triangles_),
    vertices(vertices_),
    areas(computeAreas(triangles_, vertices_)),
    stiffnessMatrix(computeStiffnessMatrix(triangles_, vertices_)),
    diagonal(stiffnessMatrix.diagonal())
{
    type = FunctionalType::DirichletEnergy;
    CORRADE_ASSERT(!Math::isNan(areas), "Dirichlet Energy : areas contains NaN",);

#ifndef NDEBUG
    Eigen::MatrixXd V(vertices.size(), 3);
    for (int j = 0; j < vertices.size(); ++j) {
        for (int i = 0; i < 3; ++i) {
            V(j,i) = vertices[j][i];
        }
    }
    Eigen::MatrixXi F(triangles.size(), 3);
    for (int k = 0; k < triangles.size(); ++k) {
        for (int j = 0; j < 3; ++j) {
            F(k,j) = triangles[k][j];
        }
    }
    Eigen::SparseMatrix<Double> L, M;
    igl::cotmatrix(V, F, L);
    CORRADE_ASSERT(!Math::isNan(Containers::ArrayView(stiffnessMatrix.valuePtr(), stiffnessMatrix.nonZeros())), "gradient matrix contains nan",);
    CORRADE_ASSERT(L.isApprox(-stiffnessMatrix), "Dirichlet Energy : stiffness matrix is wrong",);
#endif
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


