//
// Created by janos on 10.12.19.
//

#pragma once

#include "GraphCommon.hpp"
#include "UniqueFunction.h"
#include "Heap.h"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/Assert.h>

#include <numeric>
#include <limits>

template<class R>
class Dijkstra {
public:

    Dijkstra() = default;

    explicit Dijkstra(R const& adjacencyList, Containers::ArrayView<const T> weights) :
            m_adjacencyList(&adjacencyList),
            m_dist(adjacencyList.size(), std::numeric_limits<double>::infinity()),
            m_prev(adjacencyList.size(), -1),
            m_weights(weights) {
    }

    void setSource(int source) { m_heap.emplace(source, 0); }

    void reset() {
        std::fill(m_dist.begin(), m_dist.end(), std::numeric_limits<double>::infinity());
        std::fill(m_prev.begin(), m_prev.end(), -1);
        m_heap = decltype(m_heap)();
    }

    bool step(int& u, double& d) {
        if(m_heap.empty()) return false;

        auto const& top = m_heap.top();
        u = top.node;
        d = top.distance;
        m_heap.pop();

        if(d > m_dist[u])
            return true;

        for(auto const& [v, w]: (*m_adjacencyList)[u]){
            auto alt = detail::detach(m_weights[w]) + detail::detach(d);
            if(alt < m_dist[v]){
                m_dist[v] = alt;
                m_prev[v] = u;
                m_heap.emplace(v, alt);
            }
        }
        return true;
    }

    template<class... F>
    auto run(F&& ... cbs) {
        int u;
        double d;
        while(step(u, d)){
            if((cbs(u) || ... ))
                break;
        }
    }

    Graph::ReversedShortestPath<Dijkstra> getShortestPathReversed(const int start, const int target) const {
        auto& self = *this;
        CORRADE_INTERNAL_ASSERT(target >= 0 && start >= 0);
        return {{target, self},
                {start,  self}};
    }

private:
    friend Graph::ReversedPathIterator<Dijkstra>;



    Heap<Graph::HeapElement> m_heap;
    R const* m_adjacencyList;
    Containers::Array<double> m_dist;
    Containers::Array<int> m_prev;
    Containers::ArrayView<const T> m_weights;
};
