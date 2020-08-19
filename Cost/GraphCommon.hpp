//
// Created by janos on 03.03.20.
//

#pragma once

#include <Corrade/Utility/Assert.h>

namespace Graph {

struct HeapElement {
    HeapElement(int n, double d) : node(n), distance(d) {}

    int node;
    double distance;

    auto operator<=>(const HeapElement& other) const {
        return distance <=> other.distance;
    }
};

struct InverseDistanceCompare {
    template<class T>
    bool operator()(T const& el1, T const& el2) {
        return el1.distance > el2.distance;
    }
};

struct Edge {
    Edge(int a_, int b_) : a(a_), b(b_) {}

    int a, b;
#ifdef __cpp_impl_three_way_comparison

    auto operator<=>(const Edge&) const = default;

#else

    bool operator<(const Edge& other) const{
        return std::tie(a, b) < std::tie(other.a, other.b);
    }

    bool operator==(const Edge& other) const{
        return std::tie(a, b) == std::tie(other.a, other.b);
    }
#endif
};

template<class A>
class ReversedPathIterator {
public:
    using value_type = Edge;
    using difference_type = int;
    using reference = value_type;
    using pointer = void;

    ReversedPathIterator(const int node, const A& algo) :
            m_e(node, node),
            m_algo(algo) {
        CORRADE_INTERNAL_ASSERT(node >= 0);
        ++(*this);
    }

    value_type operator*() {
        return m_e;
    }

    ReversedPathIterator& operator++() {
        CORRADE_INTERNAL_ASSERT(m_e.a >= 0);
        m_e = Edge(m_algo->m_prev[m_e.a], m_e.a);
        return *this;
    }

    bool operator!=(const ReversedPathIterator& other) const {
        return other.m_e.a != m_e.a || other.m_e.b != m_e.b;
    }

private:
    Edge m_e;
    Corrade::Containers::Reference<std::add_const_t<A>> m_algo;
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