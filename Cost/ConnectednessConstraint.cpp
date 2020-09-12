//
// Created by janos on 21.05.20.
//

#include "Functional.hpp"
#include "ConnectednessConstraint.hpp"

#include <adolc/adouble.h>

namespace Phasefield {


template<class Scalar>
void ConnectednessConstraint::operator()(ArrayView<const Scalar> const& parameters,
                                                         ArrayView<const Scalar> const& weights,
                                                         Scalar& out,
                                                         ArrayView<Scalar> const& gradP,
                                                         ArrayView<Scalar> const& gradW) {

    ScopedTimer timer("Connectedness");

    //Mg::Double a,b, pathThickness;
    //Mg::Color3 pathColor;
    //bool updateComponents, updateWs, updateGrad, generateLineStrips;
    //{
    //    std::lock_guard l(*mutex);
    //    generateLineStrips = generateLineStrips_;
    //    updateComponents = bool (flags & VisualizationFlag::ConnectedComponents);
    //    updateWs = bool ( flags & VisualizationFlag::GeodesicWeights );
    //    updateGrad = bool ( flags & VisualizationFlag::Gradient );
    //    pathThickness = pathThickness_;
    //}

    F<Scalar> f{a, b};
    W bump(a, b);
    WGrad bumpGrad(a, b);

    std::size_t numFaces = triangles().size();
    std::size_t numVertices = vertices.size();

    Array<double> ws(NoInit, numFaces);
    Array<bool> inInterface(NoInit, numFaces);
    Array<Scalar> uT(NoInit, numFaces);
    for(std::size_t i = 0; i < numFaces; i++) {
        auto const& t = triangles()[i];
        uT[i] = 1./3.*(U[t[0]] + U[t[1]] + U[t[2]]);
        inInterface[i] = a <= detail::detach(uT[i]) && detail::detach(uT[i]) <= b;
        ws[i] = inInterface[i] ? bump(uT[i]) : -1.;
    }

    UnionFind set(numFaces);
    for(auto[dualV1, dualV2] : dualEdges) {
        auto fuT1 = f.eval(uT[dualV1]);
        auto fuT2 = f.eval(uT[dualV2]);
        auto w = .5*(diams[dualV1] + diams[dualV2])*.5*(fuT1 + fuT2);
        updateWeight(dualV1, w, adjacencyList[dualV2]);
        updateWeight(dualV2, w, adjacencyList[dualV1]);
        if(Math::abs(detail::detach(w)) < std::numeric_limits<double>::epsilon())
            set.unite(dualV1, dualV2);
    }

    arrayResize(roots, 0);
    arrayResize(components, NoInit, numFaces);
    for(std::size_t i = 0; i < components.size(); ++i) {
        if(inInterface[i]) {
            components[i] = set.find(i);
            arrayAppend(roots, components[i]);
        } else
            components[i] = -1;
    }

    std::sort(roots.begin(), roots.end());
    numComponents = std::unique(roots.begin(), roots.end()) - roots.begin();
    arrayResize(roots, numComponents);

    Mg::Debug{} << "Phase [" << a << "," << b << "] has " << numComponents << "connected components";

    if(numComponents <= 1) {
        *cost = 0.;
        if(jacobian)
            std::fill_n(*jacobian, numVertices, 0.);
    } else {

        Array<Scalar> W(DirectInit, numComponents, 0.);

        for(int i = 0; i < numFaces; ++i) {
            auto& c = components[i];
            if(c != -1) {
                auto it = std::lower_bound(roots.begin(), roots.end(), c);
                CORRADE_INTERNAL_ASSERT(it != roots.end() && *it == c);
                auto k = it - roots.begin();
                CORRADE_INTERNAL_ASSERT(std::abs(-1 - detail::detach(ws[i])) > 1e-6);
                W[k] += ws[i]*areas[i];
                c = k;
            }
        }

        Mg::Debug{} << W;

        //run dijkstra from each connected component except last one
        arrayResize(dijkstras, numComponents - 1);
        arrayResize(stops, numComponents - 1);
        {
            ScopedTimer t("dijkstra");
            for(std::size_t i = 0; i < numComponents - 1; ++i) {
                StoppingCriteria stop(roots[i], numComponents, components);
                Dijkstra dijk(adjacencyList);
                dijk.setSource(roots[i]);
                dijk.run(stop);
                CORRADE_ASSERT(stop.foundAll(), "Connectedness : dijkstra did not find all components", false);
                dijkstras[i] = std::move(dijk);
                stops[i] = std::move(stop);
            }
        }

        Scalar d = 0;
        Array<Scalar> distancesData(ValueInit, numComponents*numComponents);
        StridedArrayView2D<Scalar> distances(distancesData, {numComponents, numComponents});

        if(jacobian)
            std::fill_n(*jacobian, numVertices, Scalar{0});

        {
            ScopedTimer tDiff("connectedness diff");

            for(std::size_t i = 0; i < numComponents; ++i) {
                for(std::size_t j = i + 1; j < numComponents; ++j) {
                    Scalar dij = 0;
                    auto Wij = W[i]*W[j];

                    for(auto&&[a, b]: dijkstras[i].getShortestPathReversed(roots[i], stops[i].target(j))) {
                        auto& av = adjacencyList[a];
                        auto it = std::find_if(av.begin(), av.end(), [b = b](const auto& n) { return b == n.vertex; });
                        CORRADE_INTERNAL_ASSERT(it != av.end());
                        dij += it->weight;

                        if(jacobian) {
                            auto lineElement = .5*(diams[a] + diams[b]);
                            const double fgrada = f.grad(detail::detach(uT[a]));
                            const double fgradb = f.grad(detail::detach(uT[b]));
                            auto weightedFGrad = detail::detach(Wij)*lineElement*.5*(1./3.);
                            for(int k = 0; k < 3; ++k) {
                                jacobian[triangles()[a][k]] += weightedFGrad*fgrada;
                                jacobian[triangles()[b][k]] += weightedFGrad*fgradb;
                            }
                        }
                    }

                    d += dij*Wij;
                    if(jacobian) {
                        distances[i][j] = dij;
                        distances[j][i] = dij;
                    }
                }
            }

            if(jacobian) {
                for(int k = 0; k < numFaces; ++k) {
                    if(components[k] < 0) //not in interface
                        continue;
                    auto wgrad = bumpGrad(uT[k]);
                    if(std::abs(wgrad) < std::numeric_limits<double>::epsilon())
                        continue;
                    std::size_t i = components[k];
                    for(std::size_t j = 0; j < numComponents; ++j) {
                        if(i == j)
                            continue;
                        double weightedGrad = distances[i][j]*W[j]*wgrad*areas[k]/3.;
                        //each incident vertex has the same influence
                        for(int l = 0; l < 3; ++l)
                            jacobian[triangles()[k][l]] += weightedGrad;
                    }
                }
            }
        }

        auto scaleFactor = 2.; //scale everything by two because we only computed half the integral
        d *= scaleFactor;
        *cost = d;

        if(jacobian) {
            for(int i = 0; i < numVertices; ++i)
                jacobian[0][i] *= scaleFactor;
        }
    }
}

template Functional::Functional(DoubleWellPotential);

}