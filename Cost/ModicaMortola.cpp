//
// Created by janos on 2/23/20.
//
#include "ModicaMortola.h"
#include "C1Functions.h"
#include "Mesh.h"
#include "Functional.hpp"

#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>

#ifdef PHASEFIELD_WITH_ADOLC
#include <adolc/adouble.h>
#endif

namespace Phasefield {

DirichletEnergy::DirichletEnergy(Mesh& m) : mesh(m) {
    mesh.requireFaceInformation();
    mesh.requireGradientOperator();
}

size_t DirichletEnergy::numParameters() const {
    return mesh.vertexCount();
}

template<class Scalar>
void DirichletEnergy::operator()(ArrayView<const Scalar> parameters,
                                 ArrayView<const Scalar> weights,
                                 Scalar& out,
                                 ArrayView<Scalar> gradP,
                                 ArrayView<Scalar> gradW) {

    for(Face face : mesh.faces()) {
        Math::Vector3<Scalar> grad{0};
        Scalar weight{0};

        for(HalfEdge he : face.halfEdges()) {
            Vertex v = he.next().tip();
            grad += parameters[v.idx]*Math::Vector3<Scalar>{mesh.gradient[he]};
            weight += weights[v.idx];
        }

        weight /= Scalar{3};

        Scalar gradNormSquared = grad.dot();
        //CORRADE_INTERNAL_ASSERT(!Math::isNan(gradNormSquared));
        out += mesh.faceArea[face]*gradNormSquared*weight;

        if(gradP) {
            for(HalfEdge he : face.halfEdges()){
                Vertex v = he.next().tip();
                gradP[v.idx] += mesh.faceArea[face]*2*Math::dot(grad, Math::Vector3<Scalar>{mesh.gradient[he]})*weight;
            }
        }
        if(gradW) {
            for(HalfEdge he : face.halfEdges()) {
                Vertex v = he.next().tip();
                gradW[v.idx] += mesh.faceArea[face]*gradNormSquared/Scalar{3};
            }
        }
    }

    //Debug{} << Math::abs(result - out);
}

AreaRegularizer::AreaRegularizer(Mesh& m) : mesh(m) {}

size_t AreaRegularizer::numParameters() const { return mesh.vertexCount(); }

template<class Scalar>
void AreaRegularizer::operator()(ArrayView<const Scalar> parameters,
                                 ArrayView<const Scalar> weights,
                                 Scalar& out,
                                 ArrayView<Scalar> gradP,
                                 ArrayView<Scalar> gradW) {
    CORRADE_ASSERT(totalArea, "Total Area not initialized", );
    double invTotalArea = 1./(*totalArea);
    Scalar integral = 0;
    SmootherStep f;
    for(Vertex vertex : mesh.vertices()) {
        size_t idx = vertex.idx;
        integral += f.eval(parameters[idx])*weights[idx]*mesh.integral[vertex];
        if(gradP) {
            gradP[idx] += f.grad(parameters[idx])*weights[idx]*mesh.integral[vertex]*invTotalArea;
        }
        if(gradW) {
            gradW[idx] += f.eval(parameters[idx])*mesh.integral[vertex];
        }
    }
    out += (integral - 0.5*(*totalArea))*invTotalArea;
}

DoubleWellPotential::DoubleWellPotential(Mesh& m) : mesh(m) { mesh.requireIntegralOperator(); }

size_t DoubleWellPotential::numParameters() const { return mesh.vertexCount(); }

template<class Scalar>
void DoubleWellPotential::operator()(ArrayView<const Scalar> parameters,
                                     ArrayView<const Scalar> weights,
                                     Scalar& cost,
                                     ArrayView<Scalar> gradP,
                                     ArrayView<Scalar> gradW) {

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

HierarchicalRegularization::HierarchicalRegularization(Mesh& m) : mesh(m) { mesh.requireIntegralOperator(); }

size_t HierarchicalRegularization::numParameters() const { return mesh.vertexCount(); }

template<class Scalar>
void HierarchicalRegularization::operator()(ArrayView<const Scalar> parameters,
                                 ArrayView<const Scalar> weights,
                                 Scalar& out,
                                 ArrayView<Scalar> gradP,
                                 ArrayView<Scalar> gradW) {
    for(Vertex vertex : mesh.vertices()) {
        const size_t idx = vertex.idx;
        const Scalar sq = 0.5*parameters[idx]*parameters[idx];
        out += (1. - weights[idx])*sq*mesh.integral[vertex];
        if(gradP) {
            gradP[idx] += (1. - weights[idx])*parameters[idx]*mesh.integral[vertex];
        }
    }
}

/* explicit instantiations */

DEFINE_FUNCTIONAL_CONSTRUCTOR(DirichletEnergy)
DEFINE_FUNCTIONAL_CONSTRUCTOR(AreaRegularizer)
DEFINE_FUNCTIONAL_CONSTRUCTOR(DoubleWellPotential)
DEFINE_FUNCTIONAL_CONSTRUCTOR(HierarchicalRegularization)

DEFINE_FUNCTIONAL_OPERATOR(DirichletEnergy, double)
DEFINE_FUNCTIONAL_OPERATOR(AreaRegularizer, double)
DEFINE_FUNCTIONAL_OPERATOR(DoubleWellPotential, double)
DEFINE_FUNCTIONAL_OPERATOR(HierarchicalRegularization, double)

#ifdef PHASEFIELD_WITH_ADOLC
DEFINE_FUNCTIONAL_OPERATOR(DirichletEnergy, adouble)
DEFINE_FUNCTIONAL_OPERATOR(AreaRegularizer, adouble)
DEFINE_FUNCTIONAL_OPERATOR(DoubleWellPotential, adouble)
DEFINE_FUNCTIONAL_OPERATOR(HierarchicalRegularization, adouble)
#endif

}