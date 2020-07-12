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
void DirichletEnergy::operator()(Scalar const* parameters, Scalar* residual, SparseMatrix* jacobians) const
{
    auto n = numParameters();
    Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>> phasefield(parameters, n);
    auto halfGrad = stiffnessMatrix.cast<Scalar>() * phasefield;
    if(jacobians){
        Eigen::Map<Eigen::VectorXd> jac(jacobians->values.data(), numParameters());
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


AreaConstraints::AreaConstraints(Containers::Array<const Mg::Vector3d> const& vs, Containers::Array<const Mg::UnsignedInt> const& is, PhasefieldTree& t) :
    vertices(vs),
    indices(is),
    tree(t)
{
    updateInternalDataStructures();
}

Mg::UnsignedInt AreaConstraints::numParameters() const {
    return tree.phasefieldData.size();
}

Mg::UnsignedInt AreaConstraints::numResiduals() const {
    return tree.numLeafs * 2;
}

void AreaConstraints::updateInternalDataStructures(){
    weights = computeIntegralOperator(triangles(), vertices);
}

template<class Scalar>
void AreaConstraints::operator()(Scalar const* params, Scalar* cost) const {

    SmootherStep smoothStep;

    auto numParams = numParameters();
    auto numRes = numResiduals();
    auto size = tree.phasefieldSize;
    auto numPhasefields = tree.nodes.size();
    auto phasefieldSize = tree.phasefieldSize;
    auto& nodes = tree.nodes;

    Containers::StridedArrayView2D<const Scalar> phasefields{{params, numParams}, {numPhasefields, phasefieldSize}};

    auto visitor = YCombinator{
            [&](auto&& visitor, PhasefieldNode& node, Containers::Array<Scalar>& p) -> void {
                auto depth = node.depth;
                auto idx = node.idx;
                bool hasChildren = node.leftChild != PhasefieldNode::None;

                Containers::Array<Scalar> p1(hasChildren ? phasefieldSize : 0);
                if(hasChildren)
                    Cr::Utility::copy(p, p1);


                for (uint32_t i = 0; i < size; ++i) {
                    Scalar pos = smoothStep.eval(phasefields[node.leftChild][i]);
                    Scalar neg = smoothStep.eval(-phasefields[node.leftChild][i]);

                    if(hasChildren){
                        p[i] *= pos;
                        p1[i] *= neg;
                    } else {
                        /* at leaf nodes we also multiply by the nodal area */
                        cost[2*idx] += pos * p[i] * weights[i];
                        cost[2*idx + 1] += neg * p[i] * weights[i];
                    }
                }

                if(hasChildren){
                    visitor(nodes[node.leftChild], p);
                    visitor(nodes[node.rightChild], p1);
                }
            }
    };

    Containers::Array<Scalar> p(phasefieldSize);
    std::fill(p.begin(), p.end(), 1.);
    std::fill_n(cost, numRes, 0.);
    visitor(tree.root(), p);
}

template void DirichletEnergy::operator()(adouble const*, adouble*, SparseMatrix*) const;
template void DirichletEnergy::operator()(double const*, double*, SparseMatrix*) const;

template void AreaRegularizer::operator()(adouble const*, adouble*, SparseMatrix*, SparseMatrix*) const;
template void AreaRegularizer::operator()(double const*, double*, SparseMatrix*, SparseMatrix*) const;

template void AreaConstraints::operator()(adouble const*, adouble*) const;
template void AreaConstraints::operator()(double const*, double*) const;

template void DoubleWellPotential::operator()(adouble const*, adouble*, SparseMatrix*, SparseMatrix*) const;
template void DoubleWellPotential::operator()(double const*, double*, SparseMatrix*, SparseMatrix*) const;
