//
// Created by janos on 2/23/20.
//
#include "ModicaMortola.h"
#include "C1Functions.h"
#include "Mesh.h"
#include "Functional.hpp"

#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>


namespace Phasefield {

DirichletEnergy::DirichletEnergy(Mesh& m) : mesh(m) {
    mesh.requireFaceAreas();
    mesh.requireGradientOperator();
}

size_t DirichletEnergy::numParameters() const {
    return mesh.vertexCount();
}

template<class Scalar>
void DirichletEnergy::operator()(ArrayView<const Scalar> const& parameters,
                                 ArrayView<const Scalar> const& weights,
                                 Scalar& out,
                                 ArrayView<Scalar> const& gradP,
                                 ArrayView<Scalar> const& gradW) {

    for(Face face : mesh.faces()) {
        Vector3d grad{0};
        Scalar weight{0};

        for(HalfEdge he : face.halfEdges()) {
            Vertex v = he.next().tip();
            grad += parameters[v.idx]*mesh.gradient[he];
            weight += weights[v.idx];
        }

        weight /= Scalar{3};

        Scalar gradNormSquared = grad.dot();
        CORRADE_INTERNAL_ASSERT(!Math::isNan(gradNormSquared));
        out += mesh.faceArea[face]*gradNormSquared*weight;

        if(gradP) {
            for(HalfEdge he : face.halfEdges()){
                Vertex v = he.next().tip();
                gradP[v.idx] += mesh.faceArea[face]*2*Math::dot(grad, mesh.gradient[he])*weight;
            }
        }
        if(gradW) {
            //for(int j = 0; j < 3; ++j) {
            //    gradW[t[j]] += mesh.faceArea[face]*gradNormSquared/Scalar{3};
            //}
            for(HalfEdge he : face.halfEdges()) {
                Vertex v = he.next().tip();
                gradW[v.idx] += mesh.faceArea[face]*gradNormSquared/Scalar{3};
            }
        }
    }
}

AreaRegularizer::AreaRegularizer(Mesh& m) : mesh(m) {}

size_t AreaRegularizer::numParameters() const { return mesh.vertexCount(); }

template<class Scalar>
void AreaRegularizer::operator()(ArrayView<const Scalar> const& parameters,
                                 ArrayView<const Scalar> const& weights,
                                 Scalar& out,
                                 ArrayView<Scalar> const& gradP,
                                 ArrayView<Scalar> const& gradW) {

    Scalar integral = 0;
    SmootherStep f;
    for(Vertex vertex : mesh.vertices()) {
        size_t idx = vertex.idx;
        integral += f.eval(parameters[idx])*weights[idx]*mesh.integral[vertex];
        if(gradP) {
            gradP[idx] += f.grad(parameters[idx])*weights[idx]*mesh.integral[vertex];
        }
        if(gradW) {
            gradW[idx] += f.eval(parameters[idx])*mesh.integral[vertex];
        }
    }
    out += integral - areaRatio*(mesh.surfaceArea*0.5);
}

DoubleWellPotential::DoubleWellPotential(Mesh& m) : mesh(m) {
    mesh.requireIntegralOperator(); }

size_t DoubleWellPotential::numParameters() const { return mesh.vertexCount(); }

template<class Scalar>
void DoubleWellPotential::operator()(ArrayView<const Scalar> const& parameters,
                                     ArrayView<const Scalar> const& weights,
                                     Scalar& cost,
                                     ArrayView<Scalar> const& gradP,
                                     ArrayView<Scalar> const& gradW) {

    DoubleWell f;
    CORRADE_INTERNAL_ASSERT(parameters.size() == weights.size());
    CORRADE_INTERNAL_ASSERT(mesh.integral.size() == weights.size());
    for(Vertex vertex : mesh.vertices()) {
        size_t idx = vertex.idx;
        CORRADE_INTERNAL_ASSERT(idx < parameters.size());
        cost += f.eval(parameters[idx])*weights[idx]*mesh.integral[vertex];
        //Debug{} << cost;
        //Debug{} << idx;
        if(gradP) {
            gradP[idx] += f.grad(parameters[idx])*weights[idx]*mesh.integral[vertex];
        }
        if(gradW) {
            gradW[idx] += f.eval(parameters[idx])*mesh.integral[idx];
        }
    }
}

/* explicit instantiations */

template void DirichletEnergy::operator()(
        ArrayView<const double> const& parameters,
        ArrayView<const double> const& weights,
        double& out,
        ArrayView<double> const& gradP,
        ArrayView<double> const& gradW);

template void AreaRegularizer::operator()(
        ArrayView<const double> const& parameters,
        ArrayView<const double> const& weights,
        double& out,
        ArrayView<double> const& gradP,
        ArrayView<double> const& gradW);


template void DoubleWellPotential::operator()(
        ArrayView<const double> const& parameters,
        ArrayView<const double> const& weights,
        double& out,
        ArrayView<double> const& gradP,
        ArrayView<double> const& gradW);


template Functional::Functional(DirichletEnergy);

template Functional::Functional(AreaRegularizer);

template Functional::Functional(DoubleWellPotential);

}