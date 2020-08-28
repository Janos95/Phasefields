//
// Created by janos on 03.03.20.
//


#pragma once

#include "GraphCommon.hpp"
#include "CircularBuffer.h"
#include "Mesh.h"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/Assert.h>


namespace Phasefield {

class Bfs {
public:

    explicit Bfs(Mesh const&, int source) :
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

    Graph::ReversedShortestPath<Bfs> getShortestPathReversed(const int start, const int target) const {
        auto& self = *this;
        CORRADE_INTERNAL_ASSERT(target >= 0 && start >= 0);
        return {{target, self},
                {start,  self}};
    }

private:
    friend Graph::ReversedPathIterator<Bfs>;

    const R& m_adjacencyList;
    Cr::Containers::Array<int> m_prev;
    Cr::Containers::Array<bool> m_visited;
    Containers::CircularBuffer<int> m_q;
};

}