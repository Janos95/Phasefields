//
// Created by janos on 09.12.19.
//

#pragma once

#include "interpolation.hpp"
#include "quadrature_ref_triangle.hpp"
#include "dijkstra.hpp"
#include "utility.hpp"
#include "union_find.hpp"
#include "bfs.hpp"
#include "detach.hpp"
#include "small_array.hpp"
#include "robin_map.h"
#include <scoped_timer/scoped_timer.hpp>

#include <Magnum/Magnum.h>

#include <ceres/first_order_function.h>

#include <numeric>

namespace Mn = Magnum;
namespace Cr = Corrade;

enum class GradientFlag : Magnum::UnsignedByte {
    Analytic = 0,
    Automatic = 1
};



template<class Scalar>
struct ConnectednessConstraint : public ceres::FirstOrderFunction
{

    ConnectednessConstraint(
            Containers::ArrayView<const Vector3> const& vertices,
            Containers::ArrayView<const Vector3ui> const& faces,
            Float epsilon,
            Float a,
            Float b);

    bool Evaluate(double const* parameters, double* cost, double* jacobian) const override ;

    int NumParameters() const override;

    struct Neighbor {
        Neighbor(int v, double w) : vertex(v), weight(w) {}

        int vertex;
        Scalar weight;
    };

    struct Edge {
        Edge(const int v1_, const int v2_) : v1(std::min(v1_, v2_)), v2(std::max(v1_, v2_)) {}

        int v1, v2;

#ifdef __cpp_impl_three_way_comparison
        auto operator<=>(const Edge&) const = default;
#else

        bool operator<(const Edge& other) const{
            return std::tie(v1, v2) < std::tie(other.v1, other.v2);
        }

        bool operator==(const Edge& other) const{
            return std::tie(v1, v2) == std::tie(other.v1, other.v2);
        }
#endif
    };

    Containers::ArrayView<const Vector3ui> triangles;
    Containers::ArrayView<const Vector3> vertices;

    Containers::Array<Edge> dualEdges;
    Containers::Array<SmallArray<3, Neighbor>> adjacencyList;

    Containers::Array<Float> lineElements;
    Containers::Array<Float> areas;
    Containers::Array<Float> diams;

    Float epsilon;
    Float a, b;

    using graph_type = Dijkstra<Containers::Array<SmallArray<3, Neighbor>>>;

    struct ShortestPathTree{
        Dijkstra<graph_type> dijkstra;
        StoppingCriteria stop;
    };

    Containers::Array<ShortestPathTree> shortestPaths;
};

template<class Scalar>
ConnectednessConstraint<Scalar>::ConnectednessConstraint(
        Containers::ArrayView<const Vector3> const& vertices_,
        Containers::ArrayView<const Vector3ui> const& triangles_,
        Float epsilon_,
        Float a_,
        Float b_):
    vertices(vertices_),
    triangles(triangles_),
    epsilon(epsilon_), a(a_), b(b_) {

    //compute edges in dual graph
    tsl::robin_map<Edge, int> edges;
    edges.reserve(3 * triangles.size());
    Containers::arrayReserve(dualEdges, 3 * triangles.size());
    for (int i = 0; i < triangles.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            auto[it, inserted] = edges.try_emplace(Edge{triangles[i][j], triangles[i][(j + 1) % 3]}, i);
            if (!inserted) {
                //both faces share edge F(i,j) - F(i,j+1mod3)
                Containers::arrayAppend(dualEdges, Containers::InPlaceInit, it->second, i);
            }
        }
    }

    Containers::arrayResize(lineElements, Containers::ValueInit, dualEdges.size());
    adjacencyList.resize(triangles.size(), SmallArray<3, Neighbor>{Containers::NoInit});
    for (auto[v1, v2] : dualEdges) {
        adjacencyList[v1].emplace_back(v2, .0);
        adjacencyList[v2].emplace_back(v1, .0);
    }

    Containers::arrayResize(diams, triangles.size());
    for (int i = 0; i < triangles.size(); ++i) {
        auto t = triangles[i];

        auto x = vertices[t[0]] - vertices[t[1]];
        auto y = vertices[t[0]] - vertices[t[2]];
        auto z = vertices[t[1]] - vertices[t[2]];

        diams[i] = std::sqrt(std::max({(x-y).dot(), (y-z).dot(), (x-z).dot()}));
        areas[i] = Math::cross(x, y).length() * .5f;
    }

#ifndef NODEBUG
    BreadthFirstSearch bfs(adjacencyList, 0);
    bfs.run();
    CORRADE_INTERNAL_ASSERT(bfs.isConnected());
#endif
}

template<class Scalar>
bool ConnectednessConstraint<Scalar>::Evaluate(double const *phasefield,
                            double *cost,
                            double *jacobian) const {

    ScopedTimer timer("Connectedness");

    F weight(a, b);
    FGrad weightGrad(a, b);
    W bump(a, b);
    WGrad bumpGrad(a, b);

    int numFaces = triangles.size();
    int numVertices = vertices.size();

    Containers::Array<Scalar> U(Containers::ValueInit, numVertices);
    std::copy_n(phasefield, numVertices, U.begin());

    if constexpr(detail::IsDiffArray<Scalar>) {
        if (jacobian)
            for (auto &u : U)
                enoki::set_requires_gradient(u);
    }

    UnionFind set(numFaces);
    Containers::Array<bool> inInterface(Containers::DirectInit, numFaces, false);
    Containers::Array<Scalar> ws(Containers::DirectInit, numFaces, -1.);
    Containers::Array<Scalar> uT(Containers::DirectInit, numFaces, Scalar(0.));

    for (int i = 0; i < numFaces; i++) {
        auto t = triangles[i];
        uT[i] = 1. / 3. * (U[t[0]] + U[t[1]] + U[t[2]]);
        bool in = a <= detail::detach(uT[i]) && detail::detach(uT[i]) <= b;
        if (in) {
            ws[i] = bump(uT[i]);
            inInterface[i] = true;
        }
    }

    for (auto [dualV1, dualV2] : dualEdges) {
        auto fuT1 = weight(uT[dualV1]);
        auto fuT2 = weight(uT[dualV2]);
        auto w = .5 * (diams[dualV1] + diams[dualV2]) * .5 * (fuT1 + fuT2);
        updateWeight(dualV1, w, adjacencyList[dualV2]);
        updateWeight(dualV2, w, adjacencyList[dualV1]);
        if (std::abs(detail::detach(w)) < std::numeric_limits<double>::epsilon())
            set.unite(dualV1, dualV2);
    }

    Containers::Array<int> components(Containers::NoInit, set.size());
    Containers::Array<int> roots;
    for (std::size_t i = 0; i < components.size(); ++i) {
        if (inInterface[i]) {
            components[i] = set.find(i);
            Containers::arrayAppend(roots, components[i]);
        } else
            components[i] = -1;
    }

    std::sort(roots.begin(), roots.end());
    Containers::arrayResize(roots, std::unique(roots.begin(), roots.end()) - roots.begin());

    auto numComponents = roots.size();
    Debug{} << "Phase [" << a << ',' << b << "] has " << numComponents << "connected components";
    if (numComponents <= 1){
        *cost = 0.;
        if(jacobian)
            std::fill_n(jacobian, numVertices, 0.);
        return true;
    }

    Containers::Array<Scalar> W(Containers::DirectInit, numComponents, 0.);

    for (int i = 0; i < numFaces; ++i) {
        auto &c = components[i];
        if (c != -1) {
            auto it = std::lower_bound(roots.begin(), roots.end(), c);
            CORRADE_INTERNAL_ASSERT(it != roots.end() && *it == c);
            auto k = std::distance(roots.begin(), it);
            CORRADE_INTERNAL_ASSERT(std::abs(-1 - detail::detach(ws[i])) > 1e-6);
            W[k] += ws[i] * areas[i];
            c = k;
        }
    }

    //run dijkstra from each connected component except last one
    Containers::arrayResize(shortestPaths, numComponents - 1);
    {
        ScopedTimer t("dijkstra");
        for (std::size_t i = 0; i < numComponents - 1; ++i) {
            StoppingCriteria stop(roots[i], numComponents, components);
            Dijkstra dijk(adjacencyList);
            dijk.setSource(roots[i]);
            dijk.run({stop});
            CORRADE_INTERNAL_ASSERT(stop.foundAll());
            shortestPaths[i].dijkstra = std::move(dijk);
            shortestPaths[i].stop = std::move(stop);
        }
    }

    Scalar d = 0;
    Containers::Array<Scalar> distancesData(Containers::ValueInit, numComponents * numComponents);
    Containers::StridedArrayView2D<Scalar> distances(distancesData, {numComponents, numComponents});

    if(jacobian && !detail::IsDiffArray<Scalar>)
        std::fill_n(jacobian, numVertices, Scalar{0});

    {
        ScopedTimer tDiff("connectedness diff");
        auto& [dijkstras, stops] = shortestPaths;

        for (std::size_t i = 0; i < numComponents; ++i) {
            for (std::size_t j = i + 1; j < numComponents; ++j) {

                Scalar dij = 0;
                auto Wij = W[i] * W[j];

                for (auto&& [a, b]: dijkstras[i].getShortestPathReversed(roots[i], stops[i].target(j))) {

                    auto &av = adjacencyList[a];
                    auto it = std::find_if(av.begin(), av.end(), [b = b](const auto &n) { return b == n.vertex; });
                    CORRADE_INTERNAL_ASSERT(it != av.end());
                    dij += it->weight;

                    if (jacobian && !detail::IsDiffArray<Scalar>) {
                        auto lineElement = .5 * (diams[a] + diams[b]);
                        const double fgrada = weightGrad(detail::detach(uT[a]));
                        const double fgradb = weightGrad(detail::detach(uT[b]));
                        auto weightedFGrad = detail::detach(Wij) * lineElement * .5 * (1. / 3.);
                        for (auto v : triangles[a])
                            jacobian[v] += weightedFGrad * fgrada;
                        for (auto v : triangles[b])
                            jacobian[v] += weightedFGrad * fgradb;
                    }
                }

                d += dij * Wij;
                if (jacobian && !detail::IsDiffArray<Scalar>) {
                    distances[i][j] = dij;
                    distances[j][i] = dij;
                }
            }
        }

        if constexpr (!detail::IsDiffArray<Scalar>) {
            if (jacobian) {
                for (int k = 0; k < numFaces; ++k) {
                    if (components[k] < 0) //not in interface
                        continue;
                    auto wgrad = bumpGrad(uT[k]);
                    if (std::abs(wgrad) < std::numeric_limits<double>::epsilon())
                        continue;
                    std::size_t i = components[k];
                    for (std::size_t j = 0; j < numComponents; ++j) {
                        if (i == j)
                            continue;
                        double weightedGrad = distances(i, j) * W[j] * wgrad * areas[k] / 3.;
                        //each incident vertex has the same influence
                        for (int l = 0; l < 3; ++l)
                            jacobian[triangles[k][l]] += weightedGrad;
                    }
                }
            }
        }
    }

    auto scaleFactor = (1. / std::pow(epsilon, 2)) * 2.; // 1/(eps * eps)
    d *= scaleFactor;
    *cost = detail::detach(d);

    if (jacobian) {
        if constexpr(detail::IsDiffArray<Scalar>) {
            Scalar::simplify_graph_();
            enoki::backward(d);
            for (int i = 0; i < numVertices; ++i) {
                jacobian[i] = enoki::gradient(U[i]);
            }
        } else {
            for (int i = 0; i < numVertices; ++i)
                jacobian[i] *= scaleFactor;
        }
    }
    return true;
}

template<class Scalar>
int ConnectednessConstraint<Scalar>::NumParameters() const {
    return vertices.size();
}

