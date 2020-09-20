//
// Created by janos on 8/25/20.
//

#include "WeakYamabe.h"
#include "Mesh.h"
#include "C1Functions.h"
#include "Functional.hpp"

#include <Corrade/Containers/Array.h>

namespace Phasefield {

WeakYamabe::WeakYamabe(Mesh& m) : mesh(&m) {
    mesh->requireGaussianCurvature();
    mesh->requireIntegralOperator();
}

/**
 * ∫ ∇s(∇χ(u)⋅ϕ + χ(u)⋅∇ϕ) - K⋅χ(u)⋅ϕ
 */
template<class Scalar>
void WeakYamabe::operator()(const ArrayView<const Scalar>& parameters, const ArrayView<const Scalar>& weights, double& out,
                       const ArrayView<Scalar>& gradP, const ArrayView<Scalar>& gradW) {

    SmoothIndicatorFunction chi;

    size_t n = mesh->vertexCount();
    auto phasefield = parameters.prefix(n);
    auto scalingFactor = parameters.slice(n, 2*n);

    VertexData<Scalar> difference{n};

    for(Face face : mesh->faces()) {
        Vector3d gradScaling{0};
        Vector3d gradChiOfU{0};

        for(HalfEdge he : face.halfEdges()) {
            Vertex v = he.next().tip();
            gradScaling += scalingFactor[v.idx]*mesh->gradient[he];
            gradChiOfU += chi.eval(phasefield[v.idx])*mesh->gradient[he];
        }

        Scalar scaleDotChi = Math::dot(gradScaling, gradChiOfU);

        for(HalfEdge he : face.halfEdges()) {
            Vertex v = he.next().tip();
            size_t idx = v.idx;

            Scalar chiOfU = chi.eval(phasefield[v.idx]);
            Vector3d gradBasis = mesh->gradient[he];

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

}
