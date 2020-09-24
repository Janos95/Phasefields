//
// Created by janos on 03.03.20.
//


#pragma once

#include "GraphCommon.h"
#include "CircularBuffer.h"
#include "Mesh.h"

#include <queue>
#include <Corrade/Utility/Assert.h>

namespace Phasefield {

template<class V>
class Bfs {
public:

    using VertexType = V;
    using EdgeType = std::conditional_t<std::is_same_v<V, Vertex>, Edge, DualEdge>;

    explicit Bfs() = default;

    explicit Bfs(Mesh const& mesh) : m_mesh(&mesh) {
        update();
    }

    void setSource(V source) {
        m_queue.emplace(source);
        m_dist[source] = 0;
    }

    void reset() {
        m_queue = std::queue<VertexType>();
        for(EdgeType& prev : m_shortestPaths) prev = EdgeType{Invalid, m_mesh};
        for(size_t& dist : m_dist) dist = Invalid;
    }

    void update() {
        size_t n = std::is_same_v<V, Face> ? m_mesh->faceCount() : m_mesh->vertexCount();
        arrayResize(m_dist, NoInit, n);
        arrayResize(m_shortestPaths, NoInit, n);
        reset();
    }

    auto getEdges(VertexType node) {
        if constexpr (std::is_same_v<VertexType, Vertex>) return node.edges();
        if constexpr (std::is_same_v<VertexType, Face>) { return node.dualEdges(); }
    }

    auto getNeighbor(VertexType v, EdgeType edge) {
        if constexpr (std::is_same_v<EdgeType, Edge>) return edge.otherVertex(v);
        if constexpr (std::is_same_v<EdgeType , DualEdge>) return edge.otherFace(v);
    }

    bool step(V& v) {
        if(m_queue.empty()) return false;

        v = m_queue.front();
        m_queue.pop();
        size_t distance = m_dist[v];

        for(EdgeType e : getEdges(v)) {
            if constexpr (std::is_same_v<EdgeType, DualEdge>) CORRADE_ASSERT(e.isValid(), "Dual Edge not valid", false);
            VertexType neighbor = getNeighbor(v, e);
            size_t relaxedDist = 1 + distance;
            if(relaxedDist < m_dist[neighbor]) {
                m_dist[neighbor] = relaxedDist;
                m_shortestPaths[neighbor] = e;
                m_queue.emplace(neighbor);
            }
        }

        return true;
    }

    template<class... F>
    inline void run(F&& ... cbs) {
        V v;
        double d;
        while(step(v, d)) {
            if((cbs(v) || ... ))
                break;
        }
    }

    Graph::ReversedShortestPath<Bfs> getShortestPathReversed(V start, V target) const {
        return {{target, this}, {start,  this}};
    }

    bool visited(VertexType v) { return m_dist[v] != Invalid; }

private:
    friend Graph::ReversedPathIterator<Bfs>;

    Mesh const* m_mesh = nullptr;

    std::queue<VertexType> m_queue;
    MeshData<V, size_t> m_dist;
    MeshData<VertexType, EdgeType> m_shortestPaths;
};

}