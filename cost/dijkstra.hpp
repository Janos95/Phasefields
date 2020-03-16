//
// Created by janos on 10.12.19.
//

#pragma once

#include "graph_common.hpp"
#include "detach.hpp"

#include <enoki/autodiff.h>

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

    explicit Dijkstra(const R& adjacencyList):
        m_adjacencyList(adjacencyList),
        m_dist(adjacencyList.size(), std::numeric_limits<double>::infinity()),
        m_prev(adjacencyList.size(), -1)
    {
    }

    auto run(int source, std::initializer_list<folly::FunctionRef<bool(int)>> cbs)
    {
        auto comp = [](const auto& el1, const auto& el2) { return el1.distance > el2.distance; };
        std::priority_queue heap(comp, std::vector{HeapElement{source, 0}});

        while(!heap.empty())
        {
            auto [u, d] = heap.top();
            heap.pop();

           if(d > m_dist[u])
               continue;

           if(std::any_of(cbs.begin(), cbs.end(),[u=u](auto& cb){ return cb(u); }))
               return;

           for(const auto& [v, w]: m_adjacencyList[u]){
                auto alt = detail::detach(w) + detail::detach(d);
                if(alt < m_dist[v]){
                   m_dist[v] = alt;
                   m_prev[v] = u;
                   heap.emplace(v, alt);
                }
            }
        }
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

    const R& m_adjacencyList;
    std::vector<double> m_dist;
    std::vector<int> m_prev;
};

