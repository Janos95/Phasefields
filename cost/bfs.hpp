//
// Created by janos on 03.03.20.
//


#pragma once

#include "graph_common.hpp"

#include <Corrade/Containers/Reference.h>
#include <Corrade/Utility/Assert.h>

#include <deque>
#include <numeric>
#include <limits>

template<class R>
class BreadthFirstSearch
{
public:

    explicit BreadthFirstSearch(const R& adjacencyList, int source):
            m_adjacencyList(adjacencyList),
            m_prev(Containers::DirectInit, adjacencyList.size(), -1),
            m_visited(Containers::DirectInit, adjacencyList.size(), false)
    {
        m_q.push_back(source);
    }

    bool step(int& v)
    {
        if(m_q.empty())
            return false;

        v = m_q.front();
        m_q.pop_front();

        for(const auto& [w, _]: m_adjacencyList[v]){
            if(!m_visited[w]){
                m_prev[w] = v;
                m_visited[w] = true;
                m_q.push_back(w);
            }
        }
        return true;
    }

    auto run()
    {
        int node;
        while(step(node))
            ;
    }

    [[nodiscard]] bool isConnected() const {
        return std::all_of(m_visited.begin(),m_visited.end(), [](const auto& x){ return x; });
    }

    graph::ReversedShortestPath<BreadthFirstSearch> getShortestPathReversed(const int start, const int target) const {
        auto &self = *this;
        CORRADE_INTERNAL_ASSERT(target >= 0 && start >= 0);
        return {{target, self},
                {start,  self}};
    }

private:
    friend graph::ReversedPathIterator<BreadthFirstSearch>;

    const R& m_adjacencyList;
    Containers::Array<int> m_prev;
    Containers::Array<bool> m_visited;
    std::deque<int> m_q;
};

