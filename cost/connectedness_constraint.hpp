//
// Created by janos on 09.12.19.
//

#pragma once

#include "dijkstra.hpp"
#include "stopping_criteria.hpp"
#include "union_find.hpp"
#include "bfs.hpp"
#include "detach.hpp"
#include "small_array.hpp"
#include "hash.h"
#include "c1_functions.hpp"
#include "types.hpp"
#include "normalizeInto.hpp"
#include "paths.hpp"

#include <scoped_timer/scoped_timer.hpp>

#include <Corrade/Utility/Algorithms.h>
#include <Magnum/Magnum.h>
#include <Magnum/MeshTools/GenerateIndices.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <Corrade/Containers/GrowableArray.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/DebugTools/ColorMap.h>

#include <ceres/first_order_function.h>

#include <mutex>
#include <numeric>

namespace Mn = Magnum;
namespace Cr = Corrade;

template<class Scalar>
struct ConnectednessMetaData : Functional::MetaData {

    ConnectednessMetaData(
        Mg::Trade::MeshData& md,
        Cr::Containers::Array<Mg::Color3ub>& fc,
        VisualizationFlags& u,
        std::mutex& m) :
            meshData(md), faceColors(fc), update(u), mutex(m)
    {
    }

    Mg::Trade::MeshData& meshData;
    Cr::Containers::Array<Mg::Color3ub>& faceColors;
    VisualizationFlags& update;
    std::mutex& mutex;

    Mg::Double a = 0.05, b = 1.;
    Mg::Double pathThickness = 0.01;

    Cr::Containers::Array<InstanceData> instanceData; //tf from cylinder to path section
    Paths* paths = nullptr;
    bool updateInstanceData = false;
    bool generateLineStrips = false;

    /* this is called holding our 'global' mutex from the gui thread */
    void updateVis() override {
        if(updateInstanceData) {
            paths->instanceData = std::move(instanceData);
            updateInstanceData = false;
        }
    }

};

template<class Scalar>
struct ConnectednessConstraint : Functional
{
    ConnectednessConstraint(
            Cr::Containers::ArrayView<const Mg::Vector3d> const&,
            Cr::Containers::ArrayView<const Mg::Vector3ui> const&,
            Cr::Containers::Pointer<ConnectednessMetaData<Scalar>>);

    bool evaluate(double const* phasefield, double* cost, double* jacobian) const override ;


    int numParameters() const override;

    struct Neighbor {
        Neighbor(int v, double w) : vertex(v), weight(w) {}

        int vertex;
        Scalar weight;
    };

    struct Edge {
        Edge(Mg::UnsignedInt v1_, Mg::UnsignedInt v2_) : v1(std::min(v1_, v2_)), v2(std::max(v1_, v2_)) {}
        Mg::UnsignedInt v1, v2;

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

    Cr::Containers::ArrayView<const Mg::Vector3ui> triangles;
    Cr::Containers::ArrayView<const Mg::Vector3d> vertices;

    Cr::Containers::Array<Edge> dualEdges;
    mutable Cr::Containers::Array<SmallArray<3, Neighbor>> adjacencyList;

    Cr::Containers::Array<Mg::Double> lineElements;
    Cr::Containers::Array<Mg::Double> areas;
    Cr::Containers::Array<Mg::Double> diams;

    using graph_type = Cr::Containers::Array<SmallArray<3, Neighbor>>;



    decltype(auto) getMetaData() const { return dynamic_cast<ConnectednessMetaData<Scalar>&>(*metaData); }
};

template<class Scalar>
ConnectednessConstraint<Scalar>::ConnectednessConstraint(
        Cr::Containers::ArrayView<const Mg::Vector3d> const& vertices_,
        Cr::Containers::ArrayView<const Mg::Vector3ui> const& triangles_,
        Cr::Containers::Pointer<ConnectednessMetaData<Scalar>> md):
    Functional(std::move(md), FunctionalType::Connectedness),
    vertices(vertices_),
    triangles(triangles_),
    adjacencyList(Cr::Containers::DirectInit, triangles_.size(), Cr::Containers::NoInit)
{
    //compute edges in dual graph
    auto hash = [](Edge const& e) noexcept { return hash_int(reinterpret_cast<uint64_t const&>(e)); };
    std::unordered_map<Edge, int, decltype(hash)> edges(3 * triangles.size(), hash);
    Cr::Containers::arrayReserve(dualEdges, 3 * triangles.size());
    for (int i = 0; i < triangles.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            auto[it, inserted] = edges.try_emplace(Edge{triangles[i][j], triangles[i][(j + 1) % 3]}, i);
            if (!inserted) {
                //both faces share edge F(i,j) - F(i,j+1mod3)
                Cr::Containers::arrayAppend(dualEdges, Cr::Containers::InPlaceInit, it->second, i);
            }
        }
    }

    Cr::Containers::arrayResize(lineElements, Cr::Containers::ValueInit, dualEdges.size());
    for (auto[v1, v2] : dualEdges) {
        adjacencyList[v1].emplace_back(v2, .0);
        adjacencyList[v2].emplace_back(v1, .0);
    }

    Cr::Containers::arrayResize(diams, triangles.size());
    Cr::Containers::arrayResize(areas, triangles.size());
    for (int i = 0; i < triangles.size(); ++i) {
        auto t = triangles[i];

        auto x = vertices[t[0]] - vertices[t[1]];
        auto y = vertices[t[0]] - vertices[t[2]];
        auto z = vertices[t[1]] - vertices[t[2]];

        diams[i] = std::sqrt(std::max({(x-y).dot(), (y-z).dot(), (x-z).dot()}));
        areas[i] = Mg::Math::cross(x, y).length() * .5f;
    }

#ifndef NODEBUG
    BreadthFirstSearch bfs(adjacencyList, 0);
    bfs.run();
    CORRADE_ASSERT(bfs.isConnected(), "Connectedness Constraint : Dual Graph not connected",);
#endif
}

template<class Scalar>
bool ConnectednessConstraint<Scalar>::evaluate(double const *U,
                            double *cost,
                            double *jacobian) const {

    ScopedTimer timer("Connectedness");

    Mg::Double a,b, pathThickness;
    Mg::Color3 pathColor;
    bool updateComponents, updateWs, updateGrad, generateLineStrips;
    auto& meta = getMetaData();
    {
        std::lock_guard l(meta.mutex);
        a = meta.a; b = meta.b;
        generateLineStrips = meta.generateLineStrips;
        updateComponents = bool (meta.flags & VisualizationFlag::ConnectedComponents);
        updateWs = bool ( meta.flags & VisualizationFlag::GeodesicWeights );
        updateGrad = bool ( meta.flags & VisualizationFlag::Gradient );
        pathThickness = meta.pathThickness;
    }

    F<Scalar> f{a, b};
    W bump(a, b);
    WGrad bumpGrad(a, b);

    int numFaces = triangles.size();
    int numVertices = vertices.size();

    Cr::Containers::Array<Mg::Double> ws(Cr::Containers::NoInit, numFaces);
    Cr::Containers::Array<bool> inInterface(Cr::Containers::NoInit, numFaces);
    Cr::Containers::Array<Scalar> uT(Cr::Containers::NoInit, numFaces);
    for (int i = 0; i < numFaces; i++) {
        auto const& t = triangles[i];
        uT[i] = 1. / 3. * (U[t[0]] + U[t[1]] + U[t[2]]);
        inInterface[i] = a <= detail::detach(uT[i]) && detail::detach(uT[i]) <= b;
        ws[i] = inInterface[i] ? bump(uT[i]) : -1.;
    }

    UnionFind set(numFaces);
    for (auto [dualV1, dualV2] : dualEdges) {
        auto fuT1 = f.eval(uT[dualV1]);
        auto fuT2 = f.eval(uT[dualV2]);
        auto w = .5 * (diams[dualV1] + diams[dualV2]) * .5 * (fuT1 + fuT2);
        updateWeight(dualV1, w, adjacencyList[dualV2]);
        updateWeight(dualV2, w, adjacencyList[dualV1]);
        if (std::abs(detail::detach(w)) < std::numeric_limits<double>::epsilon())
            set.unite(dualV1, dualV2);
    }

    Cr::Containers::Array<int> roots;
    Cr::Containers::Array<int> components(Cr::Containers::NoInit, numFaces);
    for (std::size_t i = 0; i < components.size(); ++i) {
        if (inInterface[i]) {
            components[i] = set.find(i);
            Cr::Containers::arrayAppend(roots, components[i]);
        } else
            components[i] = -1;
    }

    std::sort(roots.begin(), roots.end());
    auto numComponents = std::unique(roots.begin(), roots.end()) - roots.begin();
    Cr::Containers::arrayResize(roots, numComponents);

    Mg::Debug{} << "Phase [" << a << "," << b << "] has " << numComponents << "connected components";

    Cr::Containers::Array<InstanceData> instanceData;
    if (numComponents <= 1){
        *cost = 0.;
        if(jacobian)
            std::fill_n(jacobian, numVertices, 0.);
    } else {

        Cr::Containers::Array<Scalar> W(Cr::Containers::DirectInit, numComponents, 0.);

        for (int i = 0; i < numFaces; ++i) {
            auto &c = components[i];
            if (c != -1) {
                auto it = std::lower_bound(roots.begin(), roots.end(), c);
                CORRADE_INTERNAL_ASSERT(it != roots.end() && *it == c);
                auto k = it - roots.begin();
                CORRADE_INTERNAL_ASSERT(std::abs(-1 - detail::detach(ws[i])) > 1e-6);
                W[k] += ws[i] * areas[i];
                c = k;
            }
        }

        Mg::Debug{} << W;

        //run dijkstra from each connected component except last one
        Cr::Containers::Array<Dijkstra<graph_type>> dijkstras(numComponents - 1);
        Cr::Containers::Array<StoppingCriteria> stops(numComponents - 1);
        {
            ScopedTimer t("dijkstra");
            for (std::size_t i = 0; i < numComponents - 1; ++i) {
                StoppingCriteria stop(roots[i], numComponents, components);
                Dijkstra dijk(adjacencyList);
                dijk.setSource(roots[i]);
                dijk.run(stop);
                CORRADE_ASSERT(stop.foundAll(), "Connectedness : dijkstra did not find all components",false);
                dijkstras[i] = std::move(dijk);
                stops[i] = std::move(stop);
            }
        }

        Scalar d = 0;
        Cr::Containers::Array<Scalar> distancesData(Cr::Containers::ValueInit, numComponents * numComponents);
        Cr::Containers::StridedArrayView2D<Scalar> distances(distancesData, {numComponents, numComponents});

        if (jacobian)
            std::fill_n(jacobian, numVertices, Scalar{0});

        Mg::Deg hue = 42.0_degf;
        {
            ScopedTimer tDiff("connectedness diff");

            for (std::size_t i = 0; i < numComponents; ++i) {
                for (std::size_t j = i + 1; j < numComponents; ++j) {
                    Mg::UnsignedInt numEdgesInPath = 0;

                    Scalar dij = 0;
                    auto Wij = W[i] * W[j];

                    auto color = Mg::Color3::fromHsv({hue += 137.5_degf, 0.75f, 0.9f});

                    for (auto&&[a, b]: dijkstras[i].getShortestPathReversed(roots[i], stops[i].target(j))) {
                        if (generateLineStrips) {
                            auto getMidPoint = [this](auto &t) {
                                return 1. / 3. * (vertices[t[0]] + vertices[t[1]] + vertices[t[2]]);
                            };
                            Mg::Vector3 v1{getMidPoint(triangles[a])}, v2{getMidPoint(triangles[b])};
                            Mg::Float l = (v2 - v1).length();
                            CORRADE_ASSERT("Connectedness Constraint : edge length is zero",
                                           (l > std::numeric_limits<Mg::Float>::epsilon()), false);
                            Mg::Vector3 dir{(v2 - v1) / l};
                            Mg::Vector3 orthogonal{dir[2], dir[2], -dir[0] - dir[1]};
                            if (orthogonal.dot() < std::numeric_limits<Mg::Float>::epsilon())
                                orthogonal = Mg::Vector3{-dir[1] - dir[2], dir[0], dir[0]};
                            orthogonal = orthogonal.normalized();
                            Mg::Matrix3 rot{orthogonal, dir, Mg::Math::cross(orthogonal, dir)};
                            CORRADE_ASSERT(rot.isOrthogonal(), "Connectedness : tf for path vis not orthogonal", false);
                            /* overestimate length slightly for better visual appeal */
                            Mg::Matrix4 scaling = Mg::Matrix4::scaling({pathThickness, 1.1 * l, pathThickness});
                            Mg::Vector3 mid{.5f * (v1 + v2)};
                            auto tf = Mg::Matrix4::from(rot, mid) * scaling;
                            Cr::Containers::arrayAppend(instanceData, Cr::Containers::InPlaceInit, tf, rot, color);
                        }

                        auto &av = adjacencyList[a];
                        auto it = std::find_if(av.begin(), av.end(), [b = b](const auto &n) { return b == n.vertex; });
                        CORRADE_INTERNAL_ASSERT(it != av.end());
                        dij += it->weight;

                        if (jacobian) {
                            auto lineElement = .5 * (diams[a] + diams[b]);
                            const double fgrada = f.grad(detail::detach(uT[a]));
                            const double fgradb = f.grad(detail::detach(uT[b]));
                            auto weightedFGrad = detail::detach(Wij) * lineElement * .5 * (1. / 3.);
                            for (int k = 0; k < 3; ++k) {
                                jacobian[triangles[a][k]] += weightedFGrad * fgrada;
                                jacobian[triangles[b][k]] += weightedFGrad * fgradb;
                            }
                        }
                    }

                    d += dij * Wij;
                    if (jacobian) {
                        distances[i][j] = dij;
                        distances[j][i] = dij;
                    }
                }
            }

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
                        double weightedGrad = distances[i][j] * W[j] * wgrad * areas[k] / 3.;
                        //each incident vertex has the same influence
                        for (int l = 0; l < 3; ++l)
                            jacobian[triangles[k][l]] += weightedGrad;
                    }
                }
            }
        }

        auto scaleFactor = 2.; //scale everything by two because we only computed half the integral
        d *= scaleFactor;
        *cost = d;

        if (jacobian) {
            for (int i = 0; i < numVertices; ++i)
                jacobian[i] *= scaleFactor;
        }
    }

    {
        std::lock_guard l(meta.mutex);
        if(updateComponents) {
            Mg::Deg hue = 42.0_degf;
            Cr::Containers::Array<Mg::Color3ub> randomColors(Cr::Containers::NoInit, numComponents);
            for (int i = 0; i < numComponents; ++i)
                randomColors[i] = Mg::Color3ub::fromHsv({hue += 137.5_degf, 0.75f, 0.9f});

            for (int i = 0; i < components.size(); ++i) {
                if (components[i] >= 0 && !randomColors.empty())
                    meta.faceColors[i] = randomColors[components[i]];
                else
                    meta.faceColors[i] = Mg::Color3ub(255,253,208);
            }

            meta.update = VisualizationFlag::ConnectedComponents;
        }
        if(updateWs) {
            auto colorMap = Mg::DebugTools::ColorMap::turbo();
            auto [min,max] = Mg::Math::minmax(ws);
            Mg::Double length = max - min;
            if(length < std::numeric_limits<Mg::Double>::epsilon())
                std::fill(meta.faceColors.begin(), meta.faceColors.end(), Mg::Color3ub(255,253,208));
            else{
                for (int i = 0; i < ws.size(); ++i) {
                    int idx = static_cast<int>((ws[i] - min) / length * (colorMap.size() - 1.));
                    CORRADE_ASSERT(0 <= idx && idx <= 255, "Bad index for colormap", false);
                    meta.faceColors[i] = colorMap[idx];
                }
            }
            meta.update = VisualizationFlag::GeodesicWeights;
        }
        if(updateGrad){
            auto coords = meta.meshData.mutableAttribute(Mg::Trade::MeshAttribute::TextureCoordinates);
            auto xcoords = Cr::Containers::arrayCast<2, Mg::Float>(coords).template slice<1>();
            normalizeInto({jacobian, numVertices}, xcoords);
            meta.update = VisualizationFlag::Gradient;
        }
        if(generateLineStrips){
            meta.instanceData = std::move(instanceData);
            meta.updateInstanceData = true;
        }
    }

    return true;
}

template<class Scalar>
int ConnectednessConstraint<Scalar>::numParameters() const {
    return vertices.size();
}

