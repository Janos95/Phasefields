//
// Created by janos on 03.03.20.
//

#pragma once

#include <Corrade/Utility/Assert.h>
#include "MeshElements.h"

namespace Phasefield::Graph {

template<class E>
struct HeapElement {

    E node;
    double distance;

    auto operator<=>(const HeapElement& other) const {
        return distance <=> other.distance;
    }
};

template<class A>
class ReversedPathIterator {
public:
    using EdgeType = typename A::EdgeType;
    using VertexType = typename A::VertexType;

    ReversedPathIterator(VertexType v, A const* algo) :
            m_e(algo->m_shortestPaths[v]),
            m_algo(algo) {
        CORRADE_INTERNAL_ASSERT(v);
        ++(*this);
    }

    EdgeType operator*() {
        return m_e;
    }

    ReversedPathIterator& operator++() {
        if constexpr(std::is_same_v<VertexType, Vertex>)
            m_e = m_algo->m_shortestPaths[m_e.vertex1()];
        else if constexpr(std::is_same_v<VertexType, Face>)
            m_e = m_algo->m_shortestPaths[m_e.face1()];
        return *this;
    }

    bool operator!=(ReversedPathIterator const& other) const {
        return m_e != other.m_e;
    }

private:
    EdgeType m_e;
    A const* m_algo;
};

template<class A>
class ReversedShortestPath {
public:
    using iterator = ReversedPathIterator<A>;

    ReversedShortestPath(const iterator& begin, const iterator& end) : m_begin(begin), m_end(end) {}

    auto begin() {
        return m_begin;
    }

    auto end() {
        return m_end;
    }

private:
    ReversedPathIterator<A> m_begin, m_end;
};





}
