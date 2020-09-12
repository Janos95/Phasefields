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
template<class E, class T>
class Dijkstra {
public:

    Dijkstra() = default;

    explicit Dijkstra(Mesh& mesh, EdgeData<T> const& weights) : m_mesh(mesh), m_weights(weights) {
        reset();
    }

    void setSource(E source) {
        m_heap.emplace(source, 0.);
        m_dist[source] = 0;
    }

    void reset() {
        m_heap.clear();
        for(size_t& prev : m_prev) prev = Invalid;
        for(double& dist : m_dist) dist = std::numeric_limits<double>::infinity();
    }

    void update() {
        size_t n = std::is_same_v<E, Face> ? m_mesh.faceCount() : m_mesh.vertexCount();
        arrayResize(m_dist, NoInit, n);
        arrayResize(m_prev, NoInit, n);
        reset();
    }

    bool step(E& node, double& distance) {
        if(m_heap.empty()) return false;

        auto top = m_heap.extractMin();
        node = top.node;
        distance = top.distance;

        if(distance > m_dist[node])
            return true;

        if constexpr(std::is_same_v<E, Face>) {
            HalfEdge he = node.halfEdge();
            HalfEdge it = he;
            do {
                E neighbor = it.twin().face();
                loopBody(node, it.edge(), neighbor, distance);
                it = it.next();
            } while (he != it);
        } else {
            for(HalfEdge he : node.outgoingHalfEdges()) {
                //Debug{} << he;
                E neighbor = he.tip();
                loopBody(node, he.edge(), neighbor, distance);
            }
        }

        return true;
    }

    void loopBody(E node, Edge e, E neighbor, double distance) {
        double weight = m_weights[e];
        double relaxedDist = weight + distance;
        if(relaxedDist < m_dist[neighbor]){
            m_dist[neighbor] = relaxedDist;
            m_prev[neighbor.idx] = node.idx;
            m_heap.emplace(neighbor, relaxedDist);
        }
    }

    template<class... F>
    inline auto run(F&& ... cbs) {
        E e;
        double d;
        while(step(e, d)) {
            if((cbs(e) || ... ))
                break;
        }

    }

    Graph::ReversedShortestPath<Dijkstra> getShortestPathReversed(E start, E target) const {
        auto& self = *this;
        return {{target, self}, {start,  self}};
    }

private:
    friend Graph::ReversedPathIterator<Dijkstra>;

    Mesh const& m_mesh;
    EdgeData<T> const& m_weights;

    Heap<Graph::HeapElement<E>> m_heap;
    MeshData<E, double> m_dist;
    Array<size_t> m_prev;
};

}