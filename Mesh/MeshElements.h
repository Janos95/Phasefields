//
// Created by janos on 9/5/20.
//

#pragma once

#include "Surface.h"
#include "Types.h"
#include "Range.h"

#include <compare>

namespace Phasefield {

struct HalfEdge {

    size_t idx;
    Mesh const* mesh;

    auto operator<=>(HalfEdge const& other) const { return idx <=> other.idx; }
    bool operator==(HalfEdge const& other) const = default;

    [[nodiscard]] HalfEdge next() const;

    [[nodiscard]] HalfEdge twin() const;

    [[nodiscard]] Vertex tip() const;

    [[nodiscard]] Vertex tail() const;

    [[nodiscard]] Face face() const;

    [[nodiscard]] Edge edge() const;

    [[nodiscard]] Vector3 direction() const;

    [[nodiscard]] HalfEdge nextOutgoingHalfEdge() const;

    [[nodiscard]] HalfEdge nextIncomingHalfEdge() const;

    [[nodiscard]] Corner corner() const;

    [[nodiscard]] Vector3d asVector() const;

    [[nodiscard]] bool isInterior() const;

    [[nodiscard]] bool onBoundaryLoop() const;

    [[nodiscard]] explicit operator bool() const { return idx != Invalid; }
};

struct Vertex {

    size_t idx;
    Mesh const* mesh;

    auto operator<=>(Vertex const& other) const { return idx <=> other.idx; }
    bool operator==(Vertex const& other) const = default;

    [[nodiscard]] double gaussianCurvature() const;

    [[nodiscard]] Vector3 const& normal() const;

    [[nodiscard]] Vector3 const& position() const;

    [[nodiscard]] Float scalar() const;

    [[nodiscard]] HalfEdge halfEdge() const;

    [[nodiscard]] VertexCornerRange corners() const;

    [[nodiscard]] VertexAdjacentVertexRange adjacentVertices() const;

    [[nodiscard]] IncomingHalfEdgeRange incomingHalfEdges() const;

    [[nodiscard]] OutgoingHalfEdgeRange outgoingHalfEdges() const;

    [[nodiscard]] explicit operator bool() const { return idx != Invalid; }

    //IncidentEdgeRange incidentEdges() const;
    //IncidentFaceRange incidentFaces() const;

};

struct Face {

    size_t idx;
    Mesh const* mesh;

    auto operator<=>(Face const& other) const { return idx <=> other.idx; }
    bool operator==(Face const& other) const = default;

    //FaceIncidentVertexRange vertices();
    //FaceIncidentFaceRange incidentFaces() const;

    [[nodiscard]] HalfEdge halfEdge() const;

    [[nodiscard]] FaceEdgeRange edges() const;

    [[nodiscard]] FaceCornerRange corners() const;

    [[nodiscard]] FaceVertexRange vertices() const;

    [[nodiscard]] FaceHalfEdgeRange halfEdges() const;

    [[nodiscard]] FaceDualEdgeRange dualEdges() const;

    [[nodiscard]] double area() const;

    [[nodiscard]] Vector3d normal() const;

    [[nodiscard]] double diameter() const;

    [[nodiscard]] explicit operator bool() const { return idx != Invalid; }

    [[nodiscard]] bool isValid() const { return idx != Invalid; }
};

struct Edge {

    size_t idx;
    Mesh const* mesh;

    auto operator<=>(Edge const& other) const { return idx <=> other.idx; }
    bool operator==(Edge const& other) const = default;

    //bool isManifold();
    [[nodiscard]] Vertex vertex1() const;

    [[nodiscard]] Vertex vertex2() const;

    [[nodiscard]] Vertex otherVertex(Vertex) const;

    [[nodiscard]] HalfEdge halfEdge() const;

    [[nodiscard]] bool hasDualEdge() const;

    [[nodiscard]] DualEdge dualEdge() const;

    [[nodiscard]] explicit operator bool() const { return idx != Invalid; }
};

struct Corner {

    size_t idx;
    Mesh const* mesh;

    auto operator<=>(Corner const& other) const { return idx <=> other.idx; }
    bool operator==(Corner const& other) const = default;

    [[nodiscard]] HalfEdge side1() const;

    [[nodiscard]] HalfEdge side2() const;

    [[nodiscard]] Vertex vertex() const;

    [[nodiscard]] Face face() const;

    [[nodiscard]] Rad angle() const;

    [[nodiscard]] HalfEdge halfEdge() const;

    [[nodiscard]] explicit operator bool() const { return idx != Invalid; }

};

struct DualEdge {
    size_t idx;
    Mesh const* mesh;

    auto operator<=>(DualEdge const& other) const { return idx <=> other.idx; }
    bool operator==(DualEdge const& other) const = default;

    [[nodiscard]] Face face1() const;

    [[nodiscard]] Face face2() const;

    [[nodiscard]] Edge edge() const;

    [[nodiscard]] Face otherFace(Face) const;

    [[nodiscard]] explicit operator bool() const { return idx != Invalid; }
};

Debug& operator<<(Debug& debug, Vertex const& v);
Debug& operator<<(Debug& debug, Face const& f);
Debug& operator<<(Debug& debug, HalfEdge const& he);
Debug& operator<<(Debug& debug, Edge const& e);
Debug& operator<<(Debug& debug, Corner const& c);
Debug& operator<<(Debug& debug, DualEdge const& de);

template<class E>
struct FaceCirculationIterator {
    HalfEdge he;
    bool justStarted = true;

    FaceCirculationIterator& operator++() {
        if constexpr (std::is_same_v<E, DualEdge>)
            do { he = he.next(); } while(!he.edge().hasDualEdge());
        else he = he.next();
        justStarted= false;
        return *this;
    }

    bool operator !=(FaceCirculationIterator const& other) const { return justStarted || he != other.he;  }

    E operator*() {
        if constexpr (std::is_same_v<E, Edge>) return he.edge();
        if constexpr (std::is_same_v<E, HalfEdge>) return he;
        if constexpr (std::is_same_v<E, Vertex>) return he.tail();
        if constexpr (std::is_same_v<E, DualEdge>) return he.edge().dualEdge();
    }
};

struct IncomingHalfEdgeIterator {
    HalfEdge he;
    bool justStarted = true;

    IncomingHalfEdgeIterator& operator++() { he = he.nextIncomingHalfEdge(); justStarted = false; return *this;}

    bool operator !=(IncomingHalfEdgeIterator const& other) const { return justStarted || he != other.he;  }

    HalfEdge operator*() const { return he; }
};

struct OutgoingHalfEdgeIterator {
    HalfEdge he;
    bool justStarted = true;

    OutgoingHalfEdgeIterator& operator++() { he = he.nextOutgoingHalfEdge(); justStarted = false; return *this;}

    bool operator !=(OutgoingHalfEdgeIterator const& other) const { return justStarted || he != other.he; }

    HalfEdge operator*() const { return he; }
};

struct CornersOfFaceIterator {
    HalfEdge he;
    bool justStarted = true;

    CornersOfFaceIterator& operator++() { he = he.next(); justStarted = false; return *this; }

    bool operator !=(CornersOfFaceIterator const& other) const { return justStarted || he != other.he; }

    Corner operator*() const { return he.corner(); }
};


struct VertexCornerIterator {
    HalfEdge he;
    bool justStarted = true;

    VertexCornerIterator& operator++() {
        justStarted = false;
        he = he.nextIncomingHalfEdge();
        if(he.onBoundaryLoop()) he = he.nextIncomingHalfEdge();
        return *this;
    }

    bool operator !=(VertexCornerIterator const& other) const { return justStarted || he != other.he; }

    Corner operator*() const { return he.corner(); }
};

struct VertexAdjacentVertexIterator {
    HalfEdge he;
    bool justStarted = true;

    VertexAdjacentVertexIterator& operator++() { he = he.nextIncomingHalfEdge(); justStarted = false; return *this; }

    bool operator !=(VertexAdjacentVertexIterator const& other) const { return justStarted || he != other.he; }

    Vertex operator*() const { return he.tail(); }
};

template<class T>
struct ElementIterator {
    T e;

    ElementIterator& operator++() {
        if constexpr (std::is_same_v<T, Corner>)
            do { ++e.idx; } while(e.idx < e.mesh->halfEdgeCount() && !isValid());
        if constexpr (std::is_same_v<T, DualEdge>)
            do { ++e.idx; } while(e.idx < e.mesh->edgeCount() && !isValid());
        else ++e.idx;
        return *this;
    }

    bool isValid() {
        if constexpr (std::is_same_v<Corner, T>) return !e.halfEdge().onBoundaryLoop();
        if constexpr (std::is_same_v<DualEdge, T>) return e.edge().hasDualEdge();
        return true;
    }

    bool operator !=(ElementIterator const& other) const { return e != other.e; }

    T operator*() const { return e; }
};

}

