//
// Created by janos on 8/24/20.
//

#include "Dijkstra.h"

#include <Corrade/Containers/Reference.h>
#include <limits>

namespace Phasefield {

Dijkstra::Dijkstra(Mesh const& mesh, Cr::Containers::ArrayView<const double> weights, Flag flag) :
    m_mesh(&mesh),
    m_dist(Cr::Containers::DirectInit, weights.size(), std::numeric_limits<double>::infinity()),
    m_prev(Cr::Containers::DirectInit, weights.size(), -1),
    m_weights(weights),
    _flag(flag)
{
}

void Dijkstra::setSource(int source) { m_heap.emplace(source, 0); }

void Dijkstra::reset() {
    for(std::size_t i = 0; i < m_dist.size(); ++i) {
        m_dist[i] = std::numeric_limits<Mg::Double>::infinity();
        m_prev[i] = -1;
    }
    m_heap.clear();
}

bool Dijkstra::step(int& id, double& distance) {
    if(m_heap.empty()) return false;

    auto top = m_heap.extractMin();
    id = top.node;
    distance = top.distance;

    for (auto [neighbor, weight] : getAdjacentNodes(id)) {
        auto newDist = distance + weight;
        if (m_prev[neighbor] != -1 && newDist < m_dist[neighbor]){
            m_prev[neighbor] = neighbor;
            m_dist[neighbor] = newDist;
            m_heap.descreaseKey(,newDist);
        }
    }

    return true;
}

Graph::ReversedShortestPath<Dijkstra> Dijkstra::getShortestPathReversed(const int start, const int target) const {
    auto& self = *this;
    CORRADE_INTERNAL_ASSERT(target >= 0 && start >= 0);
    return {{target, self},
            {start,  self}};
}

int Dijkstra::getAdjacentNodes(int idx) {
    return 0;
}

}