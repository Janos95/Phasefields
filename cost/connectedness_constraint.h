//
// Created by janos on 21.05.20.
//

#pragma once

#include "small_array.hpp"
#include "paths.hpp"
#include "functional.hpp"
#include "dijkstra.hpp"

#include <Corrade/Containers/Array.h>

namespace Mn = Magnum;
namespace Cr = Corrade;

template<class Scalar>
struct ConnectednessConstraint final : FunctionalD
{
    ConnectednessConstraint(
            Array<const Mg::Vector3d> const&,
            Array<const Mg::UnsignedInt> const&);

    auto triangles() { return arrayCast<Mg::Vector3ui>(indices); }

    bool evaluate(double const* phasefield, double* cost, double** jacobian) const override ;
    void updateInternalDataStructures() override;
    uint32_t numParameters() const override;
    uint32_t numResiduals() const override;

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

    Array<const Mg::UnsignedInt> const& indices;
    Array<const Mg::Vector3d> const& vertices;
    Array<Edge> dualEdges;

    Array<Mg::Double> lineElements;
    Array<Mg::Double> areas;
    Array<Mg::Double> diams;

    using graph_type = Array<SmallArray<3, Neighbor>>;

    mutable Array<SmallArray<3, Neighbor>> adjacencyList;
    mutable Array<Dijkstra<graph_type>> dijkstras;
    mutable Array<StoppingCriteria> stops;
    mutable Array<int> roots;
    mutable std::size_t numComponents;
    mutable Array<int> components;


    /* this is called holding our 'global' mutex from the gui thread */
    void updateVisualization(VisualizationProxy&) override;
    void drawImGuiOptions(bool&, DrawableType&, bool&) override;

    /* collects data for the updateVisualization callback
     * Always called after a call to evaluate.
     */
    void collectVisualizationData(ArrayView<const double> const& grad);

    VisualizationFlags* update;

    Mg::Double a = 0.05, b = 1.;
    Mg::Double pathThickness = 0.01;

    Cr::Containers::Array<InstanceData> instanceData; //tf from cylinder to path section
    Paths* paths = nullptr;
    bool updateInstanceData = false;
    bool generateLineStrips = false;
    bool updateGrad = false;
    bool updateComponents = false;
    bool updateWs = false;
};


