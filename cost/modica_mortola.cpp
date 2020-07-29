//
// Created by janos on 2/23/20.
//
#include "modica_mortola.hpp"
#include "SparseMatrix.h"
#include "fem.hpp"

#include <Corrade/Utility/Assert.h>
#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/Array.h>

#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <numeric>

using namespace Magnum;
using namespace Corrade;

namespace {

auto triangles(Containers::ArrayView<const UnsignedInt> const& indices) {
    return Containers::arrayCast<const Mg::Vector3ui>(indices);
}

}

DirichletEnergy::DirichletEnergy(
        Containers::Array<const Vector3d> const& vs,
        Containers::Array<const UnsignedInt> const& is) :
        indices(is),
        vertices(vs) {
    updateInternalDataStructures();
    CORRADE_ASSERT(!Math::isNan(areas), "Dirichlet Energy : areas contains NaN",);
}


uint32_t DirichletEnergy::numParameters() const {
    return vertices.size();
}

void DirichletEnergy::updateInternalDataStructures() {

    areas = computeAreas(triangles(indices), vertices);
    edgeNormalsData = gradient(triangles(indices), vertices);
    edgeNormals = Containers::StridedArrayView2D<Mg::Vector3d>{edgeNormalsData, {triangles(indices).size(), 3}};
}

template<class Scalar>
void DirichletEnergy::operator()(Containers::ArrayView<const Scalar> const& parameters,
                                 Containers::ArrayView<const Scalar> const& weights,
                                 Scalar& out,
                                 Containers::ArrayView<Scalar> const& gradP,
                                 Containers::ArrayView<Scalar> const& gradW) {

    const auto ts = triangles(indices);

    for(std::size_t i = 0; i < ts.size(); ++i){
        auto const& t = ts[i];

        Vector3d grad{0};
        Scalar weight{0};
        for(int j = 0; j < 3; ++j){
            grad += parameters[t[j]]*edgeNormals[i][j];
            weight += weights[t[j]];
        }
        weight /= Scalar{3};

        Scalar gradNormSquared = grad.dot();
        out += areas[i]*gradNormSquared*weight;

        if(gradP){
            for(int j = 0; j < 3; ++j){
                gradP[t[j]] += areas[i]*2*Math::dot(grad, edgeNormals[i][j])*weight;
            }
        }
        if(gradW){
            for(int j = 0; j < 3; ++j){
                gradW[t[j]] += areas[i]*gradNormSquared/Scalar{3};
            }
        }
    }
}

AreaRegularizer::AreaRegularizer(
        Containers::Array<const Magnum::Vector3d> const& vs,
        Containers::Array<const Mg::UnsignedInt> const& is) : indices(is), vertices(vs) {
    updateInternalDataStructures();
}

void AreaRegularizer::updateInternalDataStructures() {
    integralOperator = computeIntegralOperator(triangles(indices), vertices);
    auto areas = computeAreas(triangles(indices), vertices);
    area = std::accumulate(areas.begin(), areas.end(), 0.);
}

template<class Scalar>
void AreaRegularizer::operator()(Containers::ArrayView<const Scalar> const& parameters,
                                 Containers::ArrayView<const Scalar> const& weights,
                                 Scalar& out,
                                 Containers::ArrayView<Scalar> const& gradP,
                                 Containers::ArrayView<Scalar> const& gradW) {

    Scalar integral = 0;
    SmootherStep f;
    for(std::size_t i = 0; i < vertices.size(); ++i){
        integral += f.eval(parameters[i])*weights[i]*integralOperator[i];
        if(gradP){
            gradP[i] += f.grad(parameters[i])*weights[i]*integralOperator[i];
        }
        if(gradW){
            gradW[i] += f.eval(parameters[i])*integralOperator[i];
        }
    }
    out += integral - areaRatio*area;
}

DoubleWellPotential::DoubleWellPotential(
        Containers::Array<const Magnum::Vector3d> const& vs,
        Containers::Array<const Mg::UnsignedInt> const& is) :
        indices(is), vertices(vs) {
}

void DoubleWellPotential::updateInternalDataStructures() {
    integralOperator = computeIntegralOperator(triangles(indices), vertices);
}

template<class Scalar>
void DoubleWellPotential::operator()(Containers::ArrayView<const Scalar> const& parameters,
                                     Containers::ArrayView<const Scalar> const& weights,
                                     Scalar& cost,
                                     Containers::ArrayView<Scalar> const& gradP,
                                     Containers::ArrayView<Scalar> const& gradW) {

    DoubleWell f;
    for(std::size_t i = 0; i < vertices.size(); ++i){
        cost += f.eval(parameters[i])*weights[i]*integralOperator[i];
        if(gradP){
            gradP[i] += f.grad(parameters[i])*weights[i]*integralOperator[i];
        }
        if(gradW){
            gradW[i] += f.eval(parameters[i])*integralOperator[i];
        }
    }
}

template void DirichletEnergy::operator()(
        Containers::ArrayView<const double> const& parameters,
        Containers::ArrayView<const double> const& weights,
        double& out,
        Containers::ArrayView<double> const& gradP,
        Containers::ArrayView<double> const& gradW);

template void AreaRegularizer::operator()(
        Containers::ArrayView<const double> const& parameters,
        Containers::ArrayView<const double> const& weights,
        double& out,
        Containers::ArrayView<double> const& gradP,
        Containers::ArrayView<double> const& gradW);


template void DoubleWellPotential::operator()(
        Containers::ArrayView<const double> const& parameters,
        Containers::ArrayView<const double> const& weights,
        double& out,
        Containers::ArrayView<double> const& gradP,
        Containers::ArrayView<double> const& gradW);

