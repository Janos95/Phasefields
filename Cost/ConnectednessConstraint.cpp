//
// Created by janos on 21.05.20.
//

#include "Functional.hpp"
#include "ConnectednessConstraint.h"
#include "C1Functions.h"
#include "UnionFind.h"
#include "StoppingCriteria.h"
#include "Mesh.h"
#include "Dijkstra.h"
#include "VisualizationProxy.h"
#include "ImGuiWidgets.h"
#include "Viewer.h"

#include <Corrade/Utility/ConfigurationGroup.h>
#include <ScopedTimer/ScopedTimer.h>

#include <imgui.h>
#include "StlAlgorithm.h"

#ifdef PHASEFIELD_WITH_ADOLC
#include <adolc/adouble.h>
#endif

namespace Phasefield {

constexpr bool verbose = true;

ConnectednessConstraint::ConnectednessConstraint(Mesh& m) : mesh(&m) {}


template<class Scalar>
void ConnectednessConstraint::operator()(ArrayView<const Scalar> parameters,
                                         ArrayView<const Scalar> weights,
                                         Scalar& cost,
                                         ArrayView<Scalar> gradP,
                                         ArrayView<Scalar> gradW) {

    ScopedTimer timer("Connectedness", verbose);

    recalculateInterval();
    double intervalBoundary = 1.5;
    if(edge0 > edge1) intervalBoundary = -1.5;
    auto [a,b] = Math::minmax(edge1, intervalBoundary);
    //W bump{a,b};
    ParametricSmoothStep bump{edge0, edge1};
    F f{a, b};

    //if constexpr (std::is_same_v<double, Scalar>) {
    //    Debug{} << edge0 << " " << edge1;
    //    Debug{} << bump.eval(edge0);
    //    Debug{} << bump.eval(edge1);
    //    Debug{} << f.eval(edge0);
    //    Debug{} << f.eval(edge1);
    //}

    WeightExitPenalty weightPenalty;

    size_t numFaces = mesh->faceCount();

    FaceData<Scalar> ws(numFaces);
    FaceData<bool> inInterface(NoInit, numFaces);
    FaceData<Scalar> uT(numFaces);
    FaceData<Scalar> facePenalties(numFaces);
    for(Face face : mesh->faces()) {

        Scalar sumU{0};
        Scalar penalty{0};
        for(Vertex v : face.vertices()) {
            sumU += parameters[v.idx];
            penalty += weights[v.idx];
        }
        uT[face] = 1./3.*sumU;
        facePenalties[face] = weightPenalty.eval<Scalar>(penalty*1./3.);

        inInterface[face] = purePhase == 1 ? uT[face] > edge1 : uT[face] < edge1;
        ws[face] = inInterface[face] ? bump.eval(uT[face]) : -1.;
    }

    UnionFind set(numFaces);

    DualEdgeData<Scalar> dualEdgeWeights{mesh->edgeCount()};

    {
        ScopedTimer t{"Computing connected components", verbose};
        for(DualEdge de : mesh->dualEdges()) {
            Face face1 = de.face1();
            Face face2 = de.face2();

            Scalar fuT1 = f.eval(uT[face1])*facePenalties[face1];
            Scalar fuT2 = f.eval(uT[face2])*facePenalties[face2];

            Scalar w = .5*(face1.diameter() + face2.diameter())*.5*(fuT1 + fuT2);
            dualEdgeWeights[de] = w;
            if(inInterface[face1] && inInterface[face2]) {
                set.unite(face1.idx, face2.idx);
            }
        }
    }

    Array<Face> roots;
    FaceData<size_t> components{NoInit, numFaces};

    for(Face face : mesh->faces()) {
        if(inInterface[face]) {
            components[face] = set.find(face.idx);
            CORRADE_ASSERT(set.find(components[face]) == set.find(components[face]), "huh",);
            arrayAppend(roots, InPlaceInit, components[face], mesh);
            CORRADE_ASSERT(inInterface[components[face]], "Set Root is not in interface",);
        } else
            components[face] = Invalid;
    }

    std::sort(roots.begin(), roots.end());
    size_t numComponents = std::unique(roots.begin(), roots.end()) - roots.begin();
    arrayResize(roots, numComponents);

    for(Face r : roots) {
        CORRADE_ASSERT(components[r] == r.idx, "huh",);
    }

    CORRADE_ASSERT(std::all_of(roots.begin(), roots.end(), [](Face f) { return !!f; }), "Connectedness : Invalid root",);

    Debug{} << "Phase" << getFormattedInterval() << "has" << numComponents << "connected components";

    if(numComponents <= 1)
        return;

    Array<Scalar> W(DirectInit, numComponents, 0.);

    {
        ScopedTimer t{"Remapping indices", verbose};
        for(Face face : mesh->faces()) {
            size_t& c = components[face];
            if(c != Invalid) {
                Face* it = std::lower_bound(roots.begin(), roots.end(), Face{c, mesh});
                CORRADE_INTERNAL_ASSERT(it != roots.end() && it->idx == c);
                size_t k = it - roots.begin();
                //CORRADE_INTERNAL_ASSERT(std::abs(-1. - ws[face]) > 1e-6);
                W[k] += ws[face]*face.area();
                c = k;
            }
        }

        if constexpr(std::is_same_v<double, Scalar>) {
            Debug{} << "W";
            for(double w : W) {
                Debug{} << w;
            }
        }
    }


    if(ignoreSmallComponents) {
        ScopedTimer t("Removing small components", verbose);

        struct ConnectedComponent {
            Face root;
            size_t idx;
            Scalar size;
        };


        Array<ConnectedComponent> ccs{numComponents};
        for(size_t i = 0; i < numComponents; ++i) {
            ccs[i] = ConnectedComponent{roots[i], i, W[i]};
            CORRADE_ASSERT(components[roots[i]] == i, "Root component not in interface", );
        }

        numComponents = 2;
        std::nth_element(ccs.begin(), ccs.begin() + numComponents, ccs.end(), [](auto& p1, auto& p2){ return p1.size > p2.size; });
        arrayResize(roots, numComponents);

        for(Face face : mesh->faces()) {
            size_t& c = components[face];
            if(c == ccs[0].idx) {
                c = 0;
            } else if(c == ccs[1].idx) {
                c = 1;
            } else {
                c = Invalid;
            }
        }

        for(size_t i = 0; i < numComponents; ++i) {
            roots[i] = ccs[i].root;
            W[i] = ccs[i].size;
            CORRADE_ASSERT(components[roots[i]] == i, "Root component not in interface", );
        }

        if(numComponents <= 1) return;

        if constexpr (std::is_same_v<Scalar, double>)
            Debug{} << W;
    }

    /* run dijkstra from each connected component except last one */
    Array<Dijkstra<Face, Scalar>> dijkstras{numComponents - 1};
    Array<StoppingCriteria> stops{numComponents - 1};
    {
        ScopedTimer t("dijkstra", verbose);
#pragma omp parallel for
        for(std::size_t i = 0; i < numComponents - 1; ++i) {
            StoppingCriteria stop(roots[i], numComponents, components);
            Dijkstra<Face, Scalar> dijk(*mesh, dualEdgeWeights);
            dijk.setSource(roots[i]);
            dijk.run(stop);
            //stop.checkIfFoundAll();
            dijkstras[i] = std::move(dijk);
            stops[i] = std::move(stop);
        }
    }

    Scalar d = 0;
    Array<Scalar> distancesData(ValueInit, numComponents*numComponents);
    StridedArrayView2D<Scalar> distances(distancesData, {numComponents, numComponents});

    {
        ScopedTimer t("connectedness running over paths", verbose);

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
                        const Scalar fgrad1 = f.grad(uT[face1]);
                        const Scalar fgrad2 = f.grad(uT[face2]);
                        auto weightedFGrad = Wij*lineElement*.5*(1./3.);

                        for(Vertex vertex : face1.vertices())
                            gradP[vertex.idx] += weightedFGrad*fgrad1*facePenalties[face1];
                        for(Vertex vertex : face2.vertices())
                            gradP[vertex.idx] += weightedFGrad*fgrad2*facePenalties[face2];
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
            ScopedTimer tDiff("connectedness diff", verbose);
            for(Face face : mesh->faces()) {
                if(components[face] == Invalid) //not in interface
                    continue;
                Scalar wgrad = bump.grad(uT[face]);
                if(wgrad*wgrad < std::numeric_limits<double>::epsilon())
                    continue;
                size_t i = components[face];
                for(size_t j = 0; j < numComponents; ++j) {
                    if(i == j)
                        continue;
                    Scalar weightedGrad = distances[i][j]*W[j]*wgrad*face.area()/3.;
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

    //Debug{} << W;

    if(gradP) {
        for(Vertex v : mesh->vertices()) gradP[v.idx] *= 2.;
    }
}

size_t ConnectednessConstraint::numParameters() const { return mesh->vertexCount(); }


void ConnectednessConstraint::drawImGuiOptions(VisualizationProxy& proxy) {
    constexpr double minS = 0.0001;
    constexpr double maxS = 3;

    ImGui::SliderScalar("s", ImGuiDataType_Double, &s, &minS, &maxS, "%.5f", ImGuiSliderFlags_Logarithmic);
    ImGui::RadioButton("Phase at -1", &purePhase, -1); ImGui::SameLine();
    ImGui::RadioButton("Phase at 1", &purePhase, 1);

    recalculateInterval(); /* just alwasy recalc, especially on startup */

    ImGui::Text("Enforcing Connected in %s", getFormattedInterval());

    if(ImGui::Checkbox("Show Connected Components", &drawComponents)) {
        if(drawComponents) {
            proxy.shaderConfig = VisualizationProxy::ShaderConfig::VertexColors;
            proxy.setCallbacks(
                    [this](Node node) { draw(node); },
                    [this] { drawComponents = false; });
        } else proxy.setDefaultCallback();
        proxy.redraw();
    }

    ImGui::SameLine();
    ImGui::Checkbox("Ignore Small Components", &ignoreSmallComponents);
}

void ConnectednessConstraint::draw(Node& node) {
    auto parameters = node.phasefield();
    //auto weights = node.temporary();

    if(drawComponents) {
        FaceData<size_t> components(mesh->faceCount());
        size_t numComponents;

        {
            ParametricSmoothStep bump{edge0, edge1};
            ParametricSmoothStep f{edge1, edge0};

            size_t faceCount = mesh->faceCount();
            FaceData<double> ws(NoInit, faceCount);
            FaceData<bool> inInterface(NoInit, faceCount);
            FaceData<double> uT(NoInit, faceCount);

            for(Face face : mesh->faces()) {

                double sumU = 0;
                for(Vertex v : face.vertices())
                    sumU += parameters[v.idx];

                uT[face] = 1./3.*sumU;
                inInterface[face] = purePhase == 1 ? uT[face] > edge1 : uT[face] < edge1;
                ws[face] = inInterface[face] ? bump.eval(uT[face]) : -1.;
            }

            UnionFind set(faceCount);

            DualEdgeData<double> dualEdgeWeights{mesh->edgeCount()};

            for(DualEdge de : mesh->dualEdges()) {
                Face face1 = de.face1();
                Face face2 = de.face2();

                double fuT1 = f.eval(uT[face1]);
                double fuT2 = f.eval(uT[face2]);

                double w = .5*(face1.diameter() + face2.diameter())*.5*(fuT1 + fuT2);
                dualEdgeWeights[de] = w;
                if(Math::abs(w) < std::numeric_limits<double>::epsilon())
                    set.unite(face1.idx, face2.idx);
            }

            Array<Face> roots;

            for(Face face : mesh->faces()) {
                if(inInterface[face]) {
                    components[face] = set.find(face.idx);
                    arrayAppend(roots, InPlaceInit, components[face], mesh);
                } else
                    components[face] = Invalid;
            }

            std::sort(roots.begin(), roots.end());
            numComponents = std::unique(roots.begin(), roots.end()) - roots.begin();
            arrayResize(roots, numComponents);

            for(Face face : mesh->faces()) {
                size_t& c = components[face];
                if(c != Invalid) {
                    Face* it = std::lower_bound(roots.begin(), roots.end(), Face{c, mesh});
                    CORRADE_INTERNAL_ASSERT(it != roots.end() && it->idx == c);
                    size_t k = it - roots.begin();
                    CORRADE_INTERNAL_ASSERT(std::abs(-1. - ws[face]) > 1e-6);
                    c = k;
                }
            }
        }
        Debug{} << "Number of connected components" << numComponents;

        auto& vertexColors = getColors(numComponents);

        for(Vertex v : mesh->vertices()) {
            float n = 0.f;
            Color4 color;
            for(Face f : v.faces()) {
                size_t c = components[f];
                if(c != Invalid) {
                    n = n + 1.f;
                    color += vertexColors[c];
                }
            }
            color *= 1.f/n;
            mesh->color(v) = color;
        }
    }
}

void ConnectednessConstraint::saveParameters(Cr::Utility::ConfigurationGroup& group) const {
    group.setValue("s", s);
    group.setValue("ignoreSmallComponents", ignoreSmallComponents);
    group.setValue("purePhase", purePhase);
}

void ConnectednessConstraint::loadParameters(Cr::Utility::ConfigurationGroup const& group) {
    s = group.value<double>("s");
    ignoreSmallComponents = group.value<bool>("ignoreSmallComponents");
    purePhase = group.value<int>("purePhase");
}

const char* ConnectednessConstraint::getFormattedInterval() {
    static char formatted[100];
    static_assert(sizeof(formatted) == 100);
    if(purePhase == -1)
        snprintf(formatted, sizeof(formatted), "(-infty, %f]", edge1);
    else
        snprintf(formatted, sizeof(formatted), "[%f, infty)", edge1);
    return formatted;
}

void ConnectednessConstraint::recalculateInterval() {
    double width = Math::pow(epsilon ? *epsilon : 0.1, s);
    if(purePhase == -1) {
        edge1 = -1 + width;
        edge0 = -1 + 2*width;
    } else {
        edge0 = 1 - 2*width;
        edge1 = 1 - width;
    }
}


DEFINE_FUNCTIONAL_CONSTRUCTOR(ConnectednessConstraint)
DEFINE_FUNCTIONAL_OPERATOR(ConnectednessConstraint, double)

#ifdef PHASEFIELD_WITH_ADOLC
DEFINE_FUNCTIONAL_OPERATOR(ConnectednessConstraint, adouble)
#endif

}