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

    ReversedPathIterator(VertexType v, A const* algo) : m_v(v), m_algo(algo) {
        CORRADE_INTERNAL_ASSERT(v);
    }

    EdgeType operator*() {
        return m_algo->m_shortestPaths[m_v];
    }

    ReversedPathIterator& operator++() {
        if constexpr(std::is_same_v<VertexType, Vertex>)
            m_v = m_algo->m_shortestPaths[m_v].otherVertex(m_v);
        else if constexpr(std::is_same_v<VertexType, Face>)
            m_v = m_algo->m_shortestPaths[m_v].otherFace(m_v);
        return *this;
    }

    bool operator!=(ReversedPathIterator const& other) const {
        return m_v != other.m_v;
    }

private:
    VertexType m_v;
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
