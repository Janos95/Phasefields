//
// Created by janos on 09.12.19.
//

#pragma once

#include "dijkstra.hpp"
#include "utility.hpp"
#include "union_find.hpp"
#include "bfs.hpp"
#include "detach.hpp"
#include "small_array.hpp"
#include "hash.h"
#include "c1_functions.hpp"
#include "types.hpp"
#include "viewer.hpp"

#include <scoped_timer/scoped_timer.hpp>

#include <Corrade/Utility/Algorithms.h>
#include <Magnum/Magnum.h>
#include <Magnum/MeshTools/GenerateIndices.h>
#include <Magnum/Math/Matrix4.h>

#include <Corrade/Containers/GrowableArray.h>

#include <ceres/first_order_function.h>

#include <numeric>
#include <ostream>
#include <sstream>


#include <assimp/Exporter.hpp>
#include <assimp/scene.h>

namespace Mn = Magnum;
namespace Cr = Corrade;

enum class GradientFlag : Magnum::UnsignedByte {
    Analytic = 0,
    Automatic = 1
};

template<class Scalar>
struct UpdateConnectednessVis;


template<class Scalar>
struct ConnectednessMetaData : Functional::MetaData {

    Viewer* viewer;
    Mg::Double a = 0.05, b = 0.95;

    Cr::Containers::Array<int> components;
    Cr::Containers::Array<Double> ws;
    Cr::Containers::Array<InstanceData> instanceData;
    Mg::Double pathThickness = 0.01;
    Mg::Color3 pathColor = Mg::Color3::green();
    int numComponents;

    solver::Status operator()(solver::IterationSummary const&) {
        std::lock_guard l(viewer->mutex);
        if(flags & VisualizationFlag::Paths){
            viewer->instanceData = std::move(instanceData);
            viewer->update |= VisualizationFlag::Paths;
        }
        if(flags & VisualizationFlag::ConnectedComponents){
            Utility::copy(components, viewer->components);
            viewer->numComponents = numComponents;
            viewer->update = VisualizationFlag::ConnectedComponents;
        }
        if(flags & VisualizationFlag::GeodesicWeights){
            Utility::copy(ws, viewer->ws);
            viewer->update = VisualizationFlag::GeodesicWeights;
        }
        return solver::Status::CONTINUE;
    }

    void reset() override {
        flags &= VisualizationFlag::Paths;
    }
};

template<class Scalar>
struct ConnectednessConstraint : Functional
{
    ConnectednessConstraint(
            Containers::ArrayView<const Vector3d> const& vertices,
            Containers::ArrayView<const Vector3ui> const& faces);

    bool evaluate(double const* parameters, double* cost, double* jacobian) const override ;

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

    Containers::ArrayView<const Vector3ui> triangles;
    Containers::ArrayView<const Vector3d> vertices;

    Containers::Array<Edge> dualEdges;
    mutable Containers::Array<SmallArray<3, Neighbor>> adjacencyList;

    Containers::Array<Mg::Double> lineElements;
    Containers::Array<Mg::Double> areas;
    Containers::Array<Mg::Double> diams;

    using graph_type = Cr::Containers::Array<SmallArray<3, Neighbor>>;

    mutable Containers::Array<bool> inInterface;
    mutable Containers::Array<Scalar> uT;

    decltype(auto) getMetaData() const { return dynamic_cast<ConnectednessMetaData<Scalar>&>(*metaData); }
};

template<class Scalar>
ConnectednessConstraint<Scalar>::ConnectednessConstraint(
        Containers::ArrayView<const Vector3d> const& vertices_,
        Containers::ArrayView<const Vector3ui> const& triangles_):
    Functional(Containers::pointer<ConnectednessMetaData<Scalar>>(), FunctionalType::Connectedness),
    vertices(vertices_),
    triangles(triangles_),
    adjacencyList(Containers::DirectInit, triangles_.size(), Containers::NoInit),
    inInterface(Containers::NoInit, triangles_.size()),
    uT(Containers::NoInit, triangles_.size())
{
    //compute edges in dual graph
    auto hash = [](Edge const& e) noexcept { return hash_int(reinterpret_cast<uint64_t const&>(e)); };
    std::unordered_map<Edge, int, decltype(hash)> edges(3 * triangles.size(), hash);
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
    for (auto[v1, v2] : dualEdges) {
        adjacencyList[v1].emplace_back(v2, .0);
        adjacencyList[v2].emplace_back(v1, .0);
    }

    Containers::arrayResize(diams, triangles.size());
    Containers::arrayResize(areas, triangles.size());
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
    CORRADE_ASSERT(bfs.isConnected(), "Connectedness Constraint : Dual Graph not connected",);
#endif
}

template<class Scalar>
bool ConnectednessConstraint<Scalar>::evaluate(double const *phasefield,
                            double *cost,
                            double *jacobian) const {

    ScopedTimer timer("Connectedness");

    Mg::Double a,b, pathThickness;
    Mg::Color3 pathColor;
    bool generateLineStrips;
    auto& meta = getMetaData();
    {
        std::lock_guard l(meta.viewer->mutex);
        a = meta.a; b = meta.b;
        generateLineStrips = static_cast<bool>(meta.flags & VisualizationFlag::Paths);
        pathColor = meta.pathColor;
        pathThickness = meta.pathThickness;
    }

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

    Containers::Array<Mg::Double> weights(Containers::NoInit, numVertices);
    for (int i = 0; i < numVertices; ++i) weights[i] = bump(U[i]);

    Containers::Array<Mg::Double> ws(Containers::NoInit, numFaces);

    for (int i = 0; i < numFaces; i++) {
        auto const& t = triangles[i];
        uT[i] = 1. / 3. * (U[t[0]] + U[t[1]] + U[t[2]]);
        inInterface[i] = a <= detail::detach(uT[i]) && detail::detach(uT[i]) <= b;
        ws[i] = inInterface[i] ? bump(uT[i]) : -1.;
    }

    UnionFind set(numFaces);
    for (auto [dualV1, dualV2] : dualEdges) {
        auto fuT1 = weight(uT[dualV1]);
        auto fuT2 = weight(uT[dualV2]);
        auto w = .5 * (diams[dualV1] + diams[dualV2]) * .5 * (fuT1 + fuT2);
        updateWeight(dualV1, w, adjacencyList[dualV2]);
        updateWeight(dualV2, w, adjacencyList[dualV1]);
        if (std::abs(detail::detach(w)) < std::numeric_limits<double>::epsilon())
            set.unite(dualV1, dualV2);
    }

    Containers::Array<int> roots;
    Containers::Array<int> components(Containers::NoInit, numFaces);
    for (std::size_t i = 0; i < components.size(); ++i) {
        if (inInterface[i]) {
            components[i] = set.find(i);
            Containers::arrayAppend(roots, components[i]);
        } else
            components[i] = -1;
    }

    std::sort(roots.begin(), roots.end());
    auto numComponents = std::unique(roots.begin(), roots.end()) - roots.begin();
    Containers::arrayResize(roots, numComponents);

    Debug{} << "Phase [" << a << "," << b << "] has " << numComponents << "connected components";
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
            auto k = it - roots.begin();
            CORRADE_INTERNAL_ASSERT(std::abs(-1 - detail::detach(ws[i])) > 1e-6);
            W[k] += ws[i] * areas[i];
            c = k;
        }
    }

    Debug{} << W;

    //run dijkstra from each connected component except last one
    Containers::Array<Dijkstra<graph_type>> dijkstras(numComponents - 1);
    Containers::Array<StoppingCriteria> stops(numComponents - 1);
    {
        ScopedTimer t("dijkstra");
        for (std::size_t i = 0; i < numComponents - 1; ++i) {
            StoppingCriteria stop(roots[i], numComponents, components);
            Dijkstra dijk(adjacencyList);
            dijk.setSource(roots[i]);
            dijk.run(stop);
            CORRADE_INTERNAL_ASSERT(stop.foundAll());
            dijkstras[i] = std::move(dijk);
            stops[i] = std::move(stop);
        }
    }

    Scalar d = 0;
    Containers::Array<Scalar> distancesData(Containers::ValueInit, numComponents * numComponents);
    Containers::StridedArrayView2D<Scalar> distances(distancesData, {numComponents, numComponents});

    if(jacobian && !detail::IsDiffArray<Scalar>)
        std::fill_n(jacobian, numVertices, Scalar{0});

    Cr::Containers::Array<InstanceData> instanceData;
    {
        ScopedTimer tDiff("connectedness diff");

        for (std::size_t i = 0; i < numComponents; ++i) {
            for (std::size_t j = i + 1; j < numComponents; ++j) {
                Mg::UnsignedInt numEdgesInPath = 0;

                Scalar dij = 0;
                auto Wij = W[i] * W[j];

                for (auto&& [a, b]: dijkstras[i].getShortestPathReversed(roots[i], stops[i].target(j))) {
                    if(generateLineStrips){
                        auto getMidPoint = [this](auto& t){ return 1. / 3. * (vertices[t[0]] + vertices[t[1]] + vertices[t[2]]); };
                        Vector3 v1{getMidPoint(triangles[a])}, v2{getMidPoint(triangles[b])};
                        Float l = (v2 - v1).length();
                        CORRADE_ASSERT("Connectedness Constraint : edge length is zero", (l > std::numeric_limits<Float>::epsilon()), false);
                        Vector3 dir{(v2 - v1) / l};
                        Vector3 orthogonal{dir[2], dir[2], -dir[0]-dir[1]};
                        if(orthogonal.dot() < std::numeric_limits<Float>::epsilon()) orthogonal = Vector3{dir[1] - dir[2], dir[0], dir[0]};
                        orthogonal = orthogonal.normalized();
                        Mg::Matrix3 rot{orthogonal, dir, Math::cross(orthogonal, dir)};
                        CORRADE_ASSERT(rot.isOrthogonal(), "Connectedness : tf for path vis not orthonormal",false);
                        Mg::Matrix3 rotScaling = {rot[0] * pathThickness, rot[1] * l, rot[2] * pathThickness};
                        Vector3 mid{.5f * (vertices[b] + vertices[a])};
                        auto tf = Mg::Matrix4::from(rotScaling, mid);
                        Containers::arrayAppend(instanceData, Containers::InPlaceInit, tf, rot, pathColor);
                    }

                    auto &av = adjacencyList[a];
                    auto it = std::find_if(av.begin(), av.end(), [b = b](const auto &n) { return b == n.vertex; });
                    CORRADE_INTERNAL_ASSERT(it != av.end());
                    dij += it->weight;

                    if (jacobian && !detail::IsDiffArray<Scalar>) {
                        auto lineElement = .5 * (diams[a] + diams[b]);
                        const double fgrada = weightGrad(detail::detach(uT[a]));
                        const double fgradb = weightGrad(detail::detach(uT[b]));
                        auto weightedFGrad = detail::detach(Wij) * lineElement * .5 * (1. / 3.);
                        for (int k = 0; k < 3; ++k) {
                            jacobian[triangles[a][k]] = weightedFGrad * fgrada;
                            jacobian[triangles[b][k]] = weightedFGrad * fgradb;
                        }
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
                        double weightedGrad = distances[i][j] * W[j] * wgrad * areas[k] / 3.;
                        //each incident vertex has the same influence
                        for (int l = 0; l < 3; ++l)
                            jacobian[triangles[k][l]] += weightedGrad;
                    }
                }
            }
        }
    }

    auto scaleFactor = 2.;
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

    if(generateLineStrips) {
        std::lock_guard l(meta.viewer->mutex);
        meta.instanceData = std::move(instanceData);
    }

    //int modifyMe = 0;
    //if(modifyMe)
    //{
    //    std::ofstream out("/tmp/phase.ply");
    //    out << "ply" << std::endl;
    //    out << "format ascii 1.0\n";
    //    out << "element vertex " << vertices.size() << '\n';
    //    out << "property float x\n";
    //    out << "property float y\n";
    //    out << "property float z\n";
    //    out << "property float u\n";
    //    out << "element face " << triangles.size() << '\n';
    //    out << "property list uchar int vertex_indices\n";
    // //   out << "property float w\n";
    // //   out << "property int c\n";
    //    out << "end_header\n";

    //    for (int i = 0; i < vertices.size(); ++i) {
    //        for (int j = 0; j < 3; ++j) {
    //            out << vertices[i][j] << ' ';
    //        }
    //        out << phasefield[i] << '\n';
    //    }

    //    for (int i = 0; i < triangles.size(); ++i) {
    //        out << "3 ";
    //        for (int j = 0; j < 3; ++j) {
    //            out << triangles[i][j] << ' ';
    //        }
    //        out /*<< ws[i] << ' ' << components[i]*/ << '\n';
    //    }
    //}

    return true;
}

template<class Scalar>
int ConnectednessConstraint<Scalar>::numParameters() const {
    return vertices.size();
}

