//
// Created by janos on 8/24/20.
//
#include "Bfs.h"

BreadthFirstSearch::BreadthFirstSearch(const R& adjacencyList, int source) :
        m_adjacencyList(adjacencyList),
        m_prev(Cr::Containers::DirectInit, adjacencyList.size(), -1),
        m_visited(Cr::Containers::DirectInit, adjacencyList.size(), false) {
    m_q.emplaceBack(source);
}

bool step(int& v) {
    if(m_q.empty())
        return false;

    v = m_q.popFront();

    for(const auto&[w, _]: m_adjacencyList[v]) {
        if(!m_visited[w]) {
            m_prev[w] = v;
            m_visited[w] = true;
            m_q.emplaceBack(w);
        }
    }
    return true;
}

void setSource(Mg::UnsignedInt source) {

}

auto run() {
    int node;
    while(step(node));
}

[[nodiscard]] bool isConnected() const {
    bool all = false;
    for(auto visited : m_visited) {
        all &= visited;
        if(!all) break;
    }
    return all;
}

Graph::ReversedShortestPath<BreadthFirstSearch> getShortestPathReversed(const int start, const int target) const {
    auto& self = *this;
    CORRADE_INTERNAL_ASSERT(target >= 0 && start >= 0);
    return {{target, self},
            {start,  self}};
}

}
