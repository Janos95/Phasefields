//
// Created by janos on 10.12.19.
//

#pragma once

#include "GraphCommon.h"
#include "Heap.h"
#include "Mesh.h"
#include "MeshElements.h"

#include <limits>
#include <Corrade/Containers/GrowableArray.h>

namespace Phasefield {

namespace Cr = Corrade;

/*
 * Dijsktra shortest path algorithm. Can either operate on the face or
 * vertex graph induced by the manifold mesh.
 */
template<class V, class T>
class Dijkstra {
public:

    using VertexType = V;
    using EdgeType = std::conditional_t<std::is_same_v<V, Vertex>, Edge, DualEdge>;

    explicit Dijkstra() = default;

    explicit Dijkstra(Mesh const& mesh, MeshData<EdgeType, T> const& weights) : m_mesh(&mesh), m_weights(&weights) {
        update();
    }

    void setSource(V source) {
        m_heap.emplace(source, 0.);
        m_dist[source] = 0;
    }

    void reset() {
        m_heap.clear();
        for(EdgeType& prev : m_shortestPaths) prev = EdgeType{Invalid, m_mesh};
        for(double& dist : m_dist) dist = std::numeric_limits<double>::infinity();
    }

    void update() {
        size_t n = std::is_same_v<V, Face> ? m_mesh->faceCount() : m_mesh->vertexCount();
        arrayResize(m_dist, NoInit, n);
        arrayResize(m_shortestPaths, NoInit, n);
        reset();
    }

    auto getEdges(VertexType node) {
        if constexpr (std::is_same_v<VertexType, Vertex>) return node.edges();
        if constexpr (std::is_same_v<VertexType, Face>) return node.dualEdges();
    }

    auto getNeighbor(VertexType v, EdgeType edge) {
        if constexpr (std::is_same_v<EdgeType, Edge>) return edge.otherVertex(v);
        if constexpr (std::is_same_v<EdgeType , DualEdge>) return edge.otherFace(v);
    }

    bool step(V& node, double& distance) {
        if(m_heap.empty()) return false;

        auto top = m_heap.extractMin();
        node = top.node;
        distance = top.distance;

        if(distance > m_dist[node])
            return true;

        for(EdgeType e : getEdges(node)) {
            VertexType neighbor = getNeighbor(node, e);
            double weight = (*m_weights)[e];
            double relaxedDist = weight + distance;
            if(relaxedDist < m_dist[neighbor]){
                m_dist[neighbor] = relaxedDist;
                m_shortestPaths[neighbor] = e;
                m_heap.emplace(neighbor, relaxedDist);
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

    Graph::ReversedShortestPath<Dijkstra> getShortestPathReversed(V start, V target) const {
        return {{target, this}, {start,  this}};
    }

private:
    friend Graph::ReversedPathIterator<Dijkstra>;

    Mesh const* m_mesh = nullptr;
    MeshData<EdgeType, T> const* m_weights = nullptr;

    Heap<Graph::HeapElement<V>> m_heap;
    MeshData<V, double> m_dist;
    MeshData<VertexType, EdgeType> m_shortestPaths;
};

}