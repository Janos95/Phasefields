//
// Created by janos on 8/25/20.
//

#include "WeakYamabe.h"
#include "Mesh.h"
#include "C1Functions.h"
#include "Functional.hpp"

#include <Corrade/Containers/Array.h>

#include <adolc/adouble.h>

namespace Phasefield {

WeakYamabe::WeakYamabe(Mesh& m) : mesh(&m) {
    mesh->requireGaussianCurvature();
    mesh->requireIntegralOperator();
}

/**
 * ∫ ∇s(∇χ(u)⋅ϕ + χ(u)⋅∇ϕ) - K⋅χ(u)⋅ϕ
 */
template<class Scalar>
void WeakYamabe::operator()(
         ArrayView<const Scalar> parameters, ArrayView<const Scalar> weights, Scalar& out,
         ArrayView<Scalar> gradP, ArrayView<Scalar> gradW) {

    SmoothIndicatorFunction chi;

    size_t n = mesh->vertexCount();
    auto phasefield = parameters.prefix(n);
    auto scalingFactor = parameters.slice(n, 2*n);

    VertexData<Scalar> difference{n};

    for(Face face : mesh->faces()) {
        Math::Vector3<Scalar> gradScaling{0};
        Math::Vector3<Scalar> gradChiOfU{0};

        for(HalfEdge he : face.halfEdges()) {
            Vertex v = he.next().tip();
            Math::Vector3<Scalar> gradBasis{mesh->gradient[he]};
            gradScaling += scalingFactor[v.idx]*gradBasis;
            gradChiOfU += chi.eval(phasefield[v.idx])*gradBasis;
        }

        Scalar scaleDotChi = Math::dot(gradScaling, gradChiOfU);

        for(HalfEdge he : face.halfEdges()) {
            Vertex v = he.next().tip();
            size_t idx = v.idx;

            Scalar chiOfU = chi.eval(phasefield[v.idx]);
            Math::Vector3<Scalar> gradBasis{mesh->gradient[he]};

            difference[v] += (scaleDotChi + chiOfU*Math::dot(gradScaling, gradBasis) - v.gaussianCurvature()*chiOfU)*1./3.*face.area();
        }
    }

    for(Scalar v : difference) {
        out += v*v;
    }
}

size_t WeakYamabe::numParameters() const { return mesh->vertexCount(); }

DEFINE_FUNCTIONAL_CONSTRUCTOR(WeakYamabe)
DEFINE_FUNCTIONAL_OPERATOR(WeakYamabe, double)
DEFINE_FUNCTIONAL_OPERATOR(WeakYamabe, adouble)

}
