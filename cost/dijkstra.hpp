//
// Created by janos on 10.12.19.
//

#pragma once

#include "graph_common.hpp"
#include "detach.hpp"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/Assert.h>

#include <folly/Function.h>

#include <queue>
#include <numeric>
#include <limits>


template<class R>
class Dijkstra
{
public:

    explicit Dijkstra(R const& adjacencyList):
        m_adjacencyList(adjacencyList),
        m_dist(adjacencyList.size(), std::numeric_limits<double>::infinity()),
        m_prev(adjacencyList.size(), -1)
    {
    }

    void setSource(int source){ m_heap.emplace(HeapElement{source, 0}); }

    void reset() {
        std::fill(m_dist.begin(), m_dist.end(), std::numeric_limits<double>::infinity());
        std::fill(m_prev.begin(), m_prev.end(), -1);
        m_heap = decltype(m_heap)();
    }

    bool step(int& u, double& d){
        if(m_heap.empty()) return false;

        auto const& top = m_heap.top();
        u = top.node;
        d = top.distance;
        m_heap.pop();

        if(d > m_dist[u])
            return true;

        for(const auto& [v, w]: m_adjacencyList[u]){
            auto alt = detail::detach(w) + detail::detach(d);
            if(alt < m_dist[v]){
                m_dist[v] = alt;
                m_prev[v] = u;
                m_heap.emplace(v, alt);
            }
        }
        return true;
    }

    auto run(std::initializer_list<folly::FunctionRef<bool(int)>> cbs)
    {
        int u;
        double d;
        while( step(u, d) && std::any_of(cbs.begin(), cbs.end(),[u](auto& cb){ return cb(u); }))
            ;
    }

    graph::ReversedShortestPath<Dijkstra> getShortestPathReversed(const int start, const int target) const {
        auto &self = *this;
        CORRADE_INTERNAL_ASSERT(target >= 0 && start >= 0);
        return {{target, self},
                {start,  self}};
    }

private:
    friend graph::ReversedPathIterator<Dijkstra>;

    struct HeapElement
    {
        HeapElement(int n, double d): node(n), distance(d){}
        int node;
        double distance;
    };

    struct InverseDistanceCompare{
        template<class T>
        bool operator()(T const& el1, T const& el2){
            return el1.distance > el2.distance;
        }
    };

    std::priority_queue<HeapElement,  std::vector<HeapElement>, InverseDistanceCompare> m_heap;
    R const& m_adjacencyList;
    std::vector<double> m_dist;
    std::vector<int> m_prev;
};

