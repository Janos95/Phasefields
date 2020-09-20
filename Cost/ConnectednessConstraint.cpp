//
// Created by janos on 21.05.20.
//

#include "Functional.hpp"
#include "ConnectednessConstraint.h"
#include "C1Functions.h"
#include "UnionFind.h"
#include "StoppingCriteria.h"

#include <ScopedTimer/ScopedTimer.h>

#include "StlAlgorithm.h"

namespace Phasefield {

ConnectednessConstraint::ConnectednessConstraint(Mesh& m) : mesh(&m) {}


template<class Scalar>
void ConnectednessConstraint::operator()(ArrayView<const Scalar> const& parameters,
                                                         ArrayView<const Scalar> const& weights,
                                                         Scalar& cost,
                                                         ArrayView<Scalar> const& gradP,
                                                         ArrayView<Scalar> const& gradW) {

    ScopedTimer timer("Connectedness");

    F<Scalar> f{a, b};
    W bump(a, b);
    WGrad bumpGrad(a, b);

    size_t numFaces = mesh->faceCount();
    size_t numVertices = mesh->vertexCount();

    FaceData<double> ws(NoInit, numFaces);
    FaceData<bool> inInterface(NoInit, numFaces);
    FaceData<Scalar> uT(NoInit, numFaces);
    for(Face face : mesh->faces()) {

        double sumU = 0;
        for(Vertex v : face.vertices())
            sumU += parameters[v.idx];

        uT[face] = 1./3.*sumU;
        inInterface[face] = a<= uT[face] && uT[face] <= b;
        ws[face] = inInterface[face] ? bump(uT[face]) : -1.;
    }

    UnionFind set(numFaces);

    DualEdgeData<Scalar> dualEdgeWeights{mesh->edgeCount()};

    for(DualEdge de : mesh->dualEdges()) {
        Face face1 = de.face1();
        Face face2 = de.face2();

        Scalar fuT1 = f.eval(uT[face1]);
        Scalar fuT2 = f.eval(uT[face2]);

        Scalar w = .5*(face1.diameter() + face2.diameter())*.5*(fuT1 + fuT2);
        dualEdgeWeights[de] = w;
        if(Math::abs(w) < std::numeric_limits<double>::epsilon())
            set.unite(face1.idx, face2.idx);
    }

    Array<Face> roots;
    FaceData<size_t> components{NoInit, numFaces};

    for(Face face : mesh->faces()) {
        if(inInterface[face]) {
            components[face] = set.find(face.idx);
            arrayAppend(roots, InPlaceInit, components[face], mesh);
        } else
            components[face] = Invalid;
    }

    std::sort(roots.begin(), roots.end());
    size_t numComponents = std::unique(roots.begin(), roots.end()) - roots.begin();
    arrayResize(roots, numComponents);

    Debug{} << "Phase [" << a << "," << b << "] has " << numComponents << "connected components";

    if(numComponents <= 1) {
        cost = 0.;
    } else {

        Array<Scalar> W(DirectInit, numComponents, 0.);

        for(Face face : mesh->faces()) {
            size_t& c = components[face];
            if(c != Invalid) {
                Face* it = std::lower_bound(roots.begin(), roots.end(), Face{c, mesh});
                CORRADE_INTERNAL_ASSERT(it != roots.end() && it->idx == c);
                size_t k = it - roots.begin();
                CORRADE_INTERNAL_ASSERT(std::abs(-1. - ws[face]) > 1e-6);
                W[k] += ws[face]*face.area();
                c = k;
            }
        }

        Mg::Debug{} << W;

        /* run dijkstra from each connected component except last one */
        Array<Dijkstra<Face, Scalar>> dijkstras{numComponents - 1};
        Array<StoppingCriteria> stops{numComponents - 1};
        {
            ScopedTimer t("dijkstra");
            for(std::size_t i = 0; i < numComponents - 1; ++i) {
                StoppingCriteria stop(roots[i], numComponents, components);
                Dijkstra<Face, Scalar> dijk(*mesh, dualEdgeWeights);
                dijk.setSource(roots[i]);
                dijk.run(stop);
                CORRADE_ASSERT(stop.foundAll(), "Connectedness : dijkstra did not find all components", );
                dijkstras[i] = std::move(dijk);
                stops[i] = std::move(stop);
            }
        }

        Scalar d = 0;
        Array<Scalar> distancesData(ValueInit, numComponents*numComponents);
        StridedArrayView2D<Scalar> distances(distancesData, {numComponents, numComponents});

        {
            ScopedTimer tDiff("connectedness diff");

            for(size_t i = 0; i < numComponents; ++i) {
                for(size_t j = i + 1; j < numComponents; ++j) {
                    Scalar dij = 0;
                    auto Wij = W[i]*W[j];

                    for(DualEdge de : dijkstras[i].getShortestPathReversed(roots[i], stops[i].target(j))) {
                        dij += dualEdgeWeights[de];
                        Face face1 = de.face1();
                        Face face2 = de.face2();

                        if(gradP) {
                            Scalar lineElement = .5*(face1.diameter() + face2.diameter());
                            const double fgrad1 = f.grad(uT[face1]);
                            const double fgrad2 = f.grad(uT[face2]);
                            auto weightedFGrad = Wij*lineElement*.5*(1./3.);

                            for(Vertex vertex : face1.vertices())
                                gradP[vertex.idx] += weightedFGrad*fgrad1;
                            for(Vertex vertex : face2.vertices())
                                gradP[vertex.idx] += weightedFGrad*fgrad2;
                        }
                    }

                    d += dij*Wij;
                    if(gradP) {
                        distances[i][j] = dij;
                        distances[j][i] = dij;
                    }
                }
            }

            if(gradP) {
                for(Face face : mesh->faces()) {
                    if(components[face] == Invalid) //not in interface
                        continue;
                    Scalar wgrad = bumpGrad(uT[face]);
                    if(Math::abs(wgrad) < std::numeric_limits<double>::epsilon())
                        continue;
                    size_t i = components[face];
                    for(size_t j = 0; j < numComponents; ++j) {
                        if(i == j)
                            continue;
                        double weightedGrad = distances[i][j]*W[j]*wgrad*face.area()/3.;
                        /* each incident vertex has the same influence */
                        for(Vertex v : face.vertices())
                            gradP[v.idx] += weightedGrad;
                    }
                }
            }
        }

        /* we scale everything by two since we computed only half the integral */
        d *= 2.;
        cost = d;

        if(gradP) {
            for(Vertex v : mesh->vertices()) gradP[v.idx] *= 2.;
        }
    }
}

size_t ConnectednessConstraint::numParameters() const { return mesh->vertexCount(); }


DEFINE_FUNCTIONAL_CONSTRUCTOR(ConnectednessConstraint)
DEFINE_FUNCTIONAL_OPERATOR(ConnectednessConstraint, double)

}