//
// Created by janos on 10.12.19.
//

#pragma once

#include "GraphCommon.hpp"
#include "Heap.h"
#include "Mesh.h"

namespace Phasefield {

namespace Cr = Corrade;

class Dijkstra {
public:

    enum class Flag : Mg::UnsignedShort {
        VertexGraph,
        DualGraph
    };

    Dijkstra() = default;

    explicit Dijkstra(Mesh const& mesh, Cr::Containers::ArrayView<const double> weights, Flag flag);

    void setSource(int source);

    void reset();

    bool step(int& u, double& d);

    template<class... F>
    inline auto run(F&& ... cbs) {
        int u;
        double d;
        while(step(u, d)) {
            if((cbs(u) || ... ))
                break;
        }
    }

    Cr::Containers::ArrayView<Graph::Ed getAdjacentNodes(int idx);

    Graph::ReversedShortestPath<Dijkstra> getShortestPathReversed(const int start, const int target) const;

private:
    friend Graph::ReversedPathIterator<Dijkstra>;

    Containers::Heap<Graph::HeapElement> m_heap;
    Mesh const* m_mesh = nullptr;
    Cr::Containers::Array<double> m_dist;
    Cr::Containers::Array<int> m_prev;
    Cr::Containers::ArrayView<const double> m_weights;

    Flag _flag;
};

}