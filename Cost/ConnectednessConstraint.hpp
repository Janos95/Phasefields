//
// Created by janos on 09.12.19.
//

#pragma once

#include "ConnectednessConstraint.h"
#include "Dijkstra.h"
#include "StoppingCriteria.h"
#include "UnionFind.h"
#include "Bfs.h"
#include "SmallArray.h"
#include "C1Functions.h"
#include "Enums.h"
#include "normalizeInto.hpp"
#include "Paths.h"
#include "ImGuiWidgets.h"
#include "VisualizationProxy.h"

#include <ScopedTimer/ScopedTimer.h>

#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/GrowableArray.h>

#include <Magnum/Magnum.h>
#include <Magnum/MeshTools/GenerateIndices.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/FunctionsBatch.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/DebugTools/ColorMap.h>
#include <Magnum/Trade/MeshData.h>

#include <ceres/first_order_function.h>

#include <imgui.h>

#include <mutex>
#include <numeric>

namespace Mn = Magnum;
namespace Cr = Corrade;


///* this is called holding our 'global' mutex from the gui thread */
//template<class Scalar>
//void ConnectednessConstraint<Scalar>::updateVisualization(VisualizationProxy& proxy) {
//    if(updateInstanceData){
//        paths->instanceData = std::move(instanceData);
//        updateInstanceData = false;
//    }
//
//    if(updateComponents){
//        Mg::Deg hue = 42.0_degf;
//        Cr::Containers::Array<Mg::Color3ub> randomColors(Cr::Containers::NoInit, numComponents);
//        for(int i = 0; i < numComponents; ++i)
//            randomColors[i] = Mg::Color3ub::fromHsv({hue += 137.5_degf, 0.75f, 0.9f});
//
//        for(int i = 0; i < components.size(); ++i){
//            if(components[i] >= 0 && !randomColors.empty())
//                (*faceColors)[i] = randomColors[components[i]];
//            else
//                (*faceColors)[i] = Mg::Color3ub(255, 253, 208);
//        }
//
//        *update = VisualizationFlag::ConnectedComponents;
//    }
//    if(updateWs){
//        auto colorMap = Mg::DebugTools::ColorMap::turbo();
//        auto[min, max] = Mg::Math::minmax(ws);
//        Mg::Double length = max - min;
//        if(length < std::numeric_limits<Mg::Double>::epsilon())
//            std::fill(faceColors->begin(), faceColors->end(), Mg::Color3ub(255, 253, 208));
//        else{
//            for(int i = 0; i < ws.size(); ++i){
//                int idx = static_cast<int>((ws[i] - min)/length*(colorMap.size() - 1.));
//                CORRADE_ASSERT(0 <= idx && idx <= 255, "Bad index for colormap", false);
//                (*faceColors)[i] = colorMap[idx];
//            }
//        }
//        *update = VisualizationFlag::GeodesicWeights;
//    }
//    if(updateGrad){
//        auto coords = meshData->mutableAttribute(Mg::Trade::MeshAttribute::TextureCoordinates);
//        auto xcoords = Cr::Containers::arrayCast<2, Mg::Float>(coords).slice<1>();
//        normalizeInto({*jacobian, numVertices}, xcoords);
//        *update = VisualizationFlag::Gradient;
//    }
//
//    if(generateLineStrips){
//        instanceData = std::move(instanceData);
//        updateInstanceData = true;
//    }
//}
//
//
///**
// * returns true if an event triggered an exclusive visualizations options.
// * Also note that we can safely read from meta.flags since the optimization thread
// * only reads from those as well (even taking the lock in case we modify)
// */
//template<class Scalar>
//void ConnectednessConstraint<Scalar>::drawImGuiOptions(bool& makeExclusive, DrawableType& type, bool& evaluateProblem) {
//    ImGui::Text(this->name.c_str());
//    dragDoubleRange2("Positive Interval", &a, &b, 0.01f, -1.f, 1.f, "Min: %.2f", "Max: %.2f", 1.f);
//
//    ImGui::BeginGroup();
//    if(ImGui::Checkbox("Visualize Shortest Paths", &paths->drawPaths)){
//        if(paths->drawPaths){
//            evaluateProblem = true;
//            generateLineStrips = true;
//        } else{
//            generateLineStrips = false;
//        }
//    }
//
//    bool visGeodesicWeights = static_cast<bool>(flags & VisualizationFlag::GeodesicWeights);
//    if(ImGui::Checkbox("Geodesic Weights", &visGeodesicWeights)){
//        if(visGeodesicWeights){
//            type = DrawableType::FaceColored;
//            makeExclusive = true;
//            flags = VisualizationFlag::GeodesicWeights;
//            evaluateProblem = true;
//        } else flags &= ~VisualizationFlag::GeodesicWeights;
//    }
//    ImGui::EndGroup();
//    ImGui::SameLine();
//    ImGui::BeginGroup();
//    bool visComponents = static_cast<bool>(flags & VisualizationFlag::ConnectedComponents);
//    if(ImGui::Checkbox("Connected Components", &visComponents)){
//        if(visComponents){
//            type = DrawableType::FaceColored;
//            makeExclusive = true;
//            flags = VisualizationFlag::ConnectedComponents;
//            evaluateProblem = true;
//        } else flags &= ~VisualizationFlag::ConnectedComponents;
//    }
//
//    bool visGradient = static_cast<bool>(flags & VisualizationFlag::Gradient);
//    if(ImGui::Checkbox("Gradient", &visGradient)){
//        if(visGradient){
//            type = DrawableType::PhongDiffuse;
//            makeExclusive = true;
//            flags = VisualizationFlag::Gradient;
//            evaluateProblem = true;
//        } else flags &= ~VisualizationFlag::Gradient;
//    }
//    ImGui::EndGroup();
//}




template<class Scalar>
void ConnectednessConstraint<Scalar>::updateInternalDataStructures() {

    arrayResize(adjacencyList, Containers::DirectInit, triangles().size(), Containers::ValueInit);

    //compute edges in dual graph
    auto hash = [](Edge const& e) noexcept { return hash_int(reinterpret_cast<uint64_t const&>(e)); };
    std::unordered_map<Edge, int, decltype(hash)> edges(3*triangles().size(), hash);

    arrayResize(dualEdges, 0);
    arrayReserve(dualEdges, 3*triangles().size());
    for(int i = 0; i < triangles().size(); ++i){
        for(int j = 0; j < 3; ++j){
            auto[it, inserted] = edges.try_emplace(Edge{triangles()[i][j], triangles()[i][(j + 1)%3]}, i);
            if(!inserted){
                //both faces share edge F(i,j) - F(i,j+1mod3)
                arrayAppend(dualEdges, Cr::Containers::InPlaceInit, it->second, i);
            }
        }
    }

    arrayResize(lineElements, Cr::Containers::ValueInit, dualEdges.size());
    for(auto[v1, v2] : dualEdges){
        adjacencyList[v1].emplace_back(v2, .0);
        adjacencyList[v2].emplace_back(v1, .0);
    }

    arrayResize(diams, triangles().size());
    arrayResize(areas, triangles().size());
    for(int i = 0; i < triangles().size(); ++i){
        auto t = triangles()[i];

        auto x = vertices[t[0]] - vertices[t[1]];
        auto y = vertices[t[0]] - vertices[t[2]];
        auto z = vertices[t[1]] - vertices[t[2]];

        diams[i] = std::sqrt(std::max({(x - y).dot(), (y - z).dot(), (x - z).dot()}));
        areas[i] = Mg::Math::cross(x, y).length()*.5f;
    }

#ifndef NODEBUG
    BreadthFirstSearch bfs(adjacencyList, 0);
    bfs.run();
    CORRADE_ASSERT(bfs.isConnected(), "Connectedness Constraint : Dual Graph not connected",);
#endif
}

template<class Scalar>
bool ConnectednessConstraint<Scalar>::evaluate(Scalar const* U,
                                               Scalar* cost,
                                               SparseMatrix* jacobian) const {

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

    Cr::Containers::Array<Mg::Double> ws(Cr::Containers::NoInit, numFaces);
    Cr::Containers::Array<bool> inInterface(Cr::Containers::NoInit, numFaces);
    Cr::Containers::Array<Scalar> uT(Cr::Containers::NoInit, numFaces);
    for(std::size_t i = 0; i < numFaces; i++){
        auto const& t = triangles()[i];
        uT[i] = 1./3.*(U[t[0]] + U[t[1]] + U[t[2]]);
        inInterface[i] = a <= detail::detach(uT[i]) && detail::detach(uT[i]) <= b;
        ws[i] = inInterface[i] ? bump(uT[i]) : -1.;
    }

    UnionFind set(numFaces);
    for(auto[dualV1, dualV2] : dualEdges){
        auto fuT1 = f.eval(uT[dualV1]);
        auto fuT2 = f.eval(uT[dualV2]);
        auto w = .5*(diams[dualV1] + diams[dualV2])*.5*(fuT1 + fuT2);
        updateWeight(dualV1, w, adjacencyList[dualV2]);
        updateWeight(dualV2, w, adjacencyList[dualV1]);
        if(std::abs(detail::detach(w)) < std::numeric_limits<double>::epsilon())
            set.unite(dualV1, dualV2);
    }

    arrayResize(roots, 0);
    arrayResize(components, Containers::NoInit, numFaces);
    for(std::size_t i = 0; i < components.size(); ++i){
        if(inInterface[i]){
            components[i] = set.find(i);
            Cr::Containers::arrayAppend(roots, components[i]);
        } else
            components[i] = -1;
    }

    std::sort(roots.begin(), roots.end());
    numComponents = std::unique(roots.begin(), roots.end()) - roots.begin();
    Cr::Containers::arrayResize(roots, numComponents);

    Mg::Debug{} << "Phase [" << a << "," << b << "] has " << numComponents << "connected components";

    if(numComponents <= 1){
        *cost = 0.;
        if(jacobian)
            std::fill_n(*jacobian, numVertices, 0.);
    } else{

        Cr::Containers::Array<Scalar> W(Cr::Containers::DirectInit, numComponents, 0.);

        for(int i = 0; i < numFaces; ++i){
            auto& c = components[i];
            if(c != -1){
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
            for(std::size_t i = 0; i < numComponents - 1; ++i){
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
        Cr::Containers::Array<Scalar> distancesData(Cr::Containers::ValueInit, numComponents*numComponents);
        Cr::Containers::StridedArrayView2D<Scalar> distances(distancesData, {numComponents, numComponents});

        if(jacobian)
            std::fill_n(*jacobian, numVertices, Scalar{0});

        {
            ScopedTimer tDiff("connectedness diff");

            for(std::size_t i = 0; i < numComponents; ++i){
                for(std::size_t j = i + 1; j < numComponents; ++j){
                    Scalar dij = 0;
                    auto Wij = W[i]*W[j];

                    for(auto&&[a, b]: dijkstras[i].getShortestPathReversed(roots[i], stops[i].target(j))){
                        auto& av = adjacencyList[a];
                        auto it = std::find_if(av.begin(), av.end(), [b = b](const auto& n) { return b == n.vertex; });
                        CORRADE_INTERNAL_ASSERT(it != av.end());
                        dij += it->weight;

                        if(jacobian){
                            auto lineElement = .5*(diams[a] + diams[b]);
                            const double fgrada = f.grad(detail::detach(uT[a]));
                            const double fgradb = f.grad(detail::detach(uT[b]));
                            auto weightedFGrad = detail::detach(Wij)*lineElement*.5*(1./3.);
                            for(int k = 0; k < 3; ++k){
                                jacobian[triangles()[a][k]] += weightedFGrad*fgrada;
                                jacobian[triangles()[b][k]] += weightedFGrad*fgradb;
                            }
                        }
                    }

                    d += dij*Wij;
                    if(jacobian){
                        distances[i][j] = dij;
                        distances[j][i] = dij;
                    }
                }
            }

            if(jacobian){
                for(int k = 0; k < numFaces; ++k){
                    if(components[k] < 0) //not in interface
                        continue;
                    auto wgrad = bumpGrad(uT[k]);
                    if(std::abs(wgrad) < std::numeric_limits<double>::epsilon())
                        continue;
                    std::size_t i = components[k];
                    for(std::size_t j = 0; j < numComponents; ++j){
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

        if(jacobian){
            for(int i = 0; i < numVertices; ++i)
                jacobian[0][i] *= scaleFactor;
        }
    }

    return true;
}


template<class Scalar>
void ConnectednessConstraint<Scalar>::collectVisualizationData(ArrayView<const double> const& grad) {
    if(generateLineStrips){
        Mg::Deg hue = 42.0_degf;
        for(std::size_t i = 0; i < numComponents; ++i){
            for(std::size_t j = i + 1; j < numComponents; ++j){

                auto color = Mg::Color3::fromHsv({hue += 137.5_degf, 0.75f, 0.9f});

                for(auto&&[a, b]: dijkstras[i].getShortestPathReversed(roots[i], stops[i].target(j))){
                    auto getMidPoint = [this](auto& t) {
                        return 1./3.*(vertices[t[0]] + vertices[t[1]] + vertices[t[2]]);
                    };
                    Mg::Vector3 v1{getMidPoint(triangles()[a])}, v2{getMidPoint(triangles()[b])};
                    Mg::Float l = (v2 - v1).length();
                    CORRADE_ASSERT("Connectedness Constraint : edge length is zero",
                                   (l > std::numeric_limits<Mg::Float>::epsilon()),);
                    Mg::Vector3 dir{(v2 - v1)/l};
                    Mg::Vector3 orthogonal{dir[2], dir[2], -dir[0] - dir[1]};
                    if(orthogonal.dot() < std::numeric_limits<Mg::Float>::epsilon())
                        orthogonal = Mg::Vector3{-dir[1] - dir[2], dir[0], dir[0]};
                    orthogonal = orthogonal.normalized();
                    Mg::Matrix3 rot{orthogonal, dir, Mg::Math::cross(orthogonal, dir)};
                    CORRADE_ASSERT(rot.isOrthogonal(), "Connectedness : tf for path vis not orthogonal",);
                    /* overestimate length slightly for better visual appeal */
                    Mg::Matrix4 scaling = Mg::Matrix4::scaling({pathThickness, 1.1*l, pathThickness});
                    Mg::Vector3 mid{.5f*(v1 + v2)};
                    auto tf = Mg::Matrix4::from(rot, mid)*scaling;
                    Cr::Containers::arrayAppend(instanceData, Cr::Containers::InPlaceInit, tf, rot, color);
                }
            }
        }
    }
}

template<class Scalar>
uint32_t ConnectednessConstraint<Scalar>::numParameters() const {
    return vertices.size();
}

template<class Scalar>
uint32_t ConnectednessConstraint<Scalar>::numResiduals() const {
    return 1;
}

