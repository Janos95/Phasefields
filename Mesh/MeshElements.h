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
    Mesh* mesh;

    auto operator<=>(HalfEdge const& other) const { return idx <=> other.idx; }

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

    [[nodiscard]] bool isInterior() const;

    [[nodiscard]] bool onBoundaryLoop() const;
};

struct Vertex {

    size_t idx;
    Mesh* mesh;

    auto operator<=>(Vertex const& other) const { return idx <=> other.idx; }

    [[nodiscard]] Float gaussianCurvature() const;

    [[nodiscard]] Vector3& normal() const;

    [[nodiscard]] Vector3& position() const;

    [[nodiscard]] Float& scalar() const;

    [[nodiscard]] HalfEdge halfEdge() const;

    [[nodiscard]] VertexCornerRange corners() const;

    [[nodiscard]] VertexAdjacentVertexRange adjacentVertices() const;

    [[nodiscard]] IncomingHalfEdgeRange incomingHalfedges() const;

    [[nodiscard]] OutgoingHalfEdgeRange outgoingHalfedges() const;

    //IncidentEdgeRange incidentEdges() const;
    //IncidentFaceRange incidentFaces() const;

};

struct Face {

    size_t idx;
    Mesh* mesh;

    auto operator<=>(Face const& other) const { return idx <=> other.idx; }

    //FaceIncidentVertexRange vertices();
    //FaceIncidentFaceRange incidentFaces() const;

    [[nodiscard]] HalfEdge halfEdge() const;

    [[nodiscard]] FaceEdgeRange edges() const;

    [[nodiscard]] CornersOfFaceRange corners() const;

    [[nodiscard]] bool isValid() const;
};

struct Edge {

    size_t idx;
    Mesh* mesh;

    auto operator<=>(Edge const& other) const { return idx <=> other.idx; }

    //bool isManifold();
    Vertex vertex1();

    Vertex vertex2();

    HalfEdge halfEdge();
};


struct Corner {

    size_t idx;
    Mesh* mesh;

    auto operator<=>(Corner const& other) const { return idx <=> other.idx; }

    [[nodiscard]] HalfEdge side1() const;

    [[nodiscard]] HalfEdge side2() const;

    [[nodiscard]] Vertex vertex() const;

    [[nodiscard]] Face face() const;

    [[nodiscard]] Rad angle() const;

    [[nodiscard]] HalfEdge halfEdge() const;

};

struct FaceEdgeIterator {
    HalfEdge he;
    bool justStarted = true;

    FaceEdgeIterator& operator++() { he = he.next(); justStarted= false; return *this;}

    bool operator !=(FaceEdgeIterator const& other) const { return justStarted || he.idx != other.he.idx;  }

    Edge operator*() const { return he.edge(); }
};

struct IncomingHalfEdgeIterator {
    HalfEdge he;
    bool justStarted = true;

    IncomingHalfEdgeIterator& operator++() { he = he.nextIncomingHalfEdge(); justStarted = false; return *this;}

    bool operator !=(IncomingHalfEdgeIterator const& other) const { return justStarted || he.idx != other.he.idx;  }

    HalfEdge operator*() const { return he; }
};

struct OutgoingHalfEdgeIterator {
    HalfEdge he;
    bool justStarted = true;

    OutgoingHalfEdgeIterator& operator++() { he = he.nextOutgoingHalfEdge(); justStarted = false; return *this;}

    bool operator !=(OutgoingHalfEdgeIterator const& other) const { return justStarted || he.idx != other.he.idx; }

    HalfEdge operator*() const { return he; }
};

struct CornersOfFaceIterator {
    HalfEdge he;
    bool justStarted = true;

    CornersOfFaceIterator& operator++() { he = he.next(); justStarted = false; return *this; }

    bool operator !=(CornersOfFaceIterator const& other) const { return justStarted || he.idx != other.he.idx; }

    Corner operator*() const { return he.corner(); }
};


struct VertexCornerIterator {
    HalfEdge he;
    bool justStarted = true;

    VertexCornerIterator& operator++() { he = he.nextIncomingHalfEdge(); justStarted = false; return *this; }

    bool operator !=(VertexCornerIterator const& other) const { return justStarted || he.idx != other.he.idx; }

    Corner operator*() const { return he.corner(); }
};

struct VertexAdjacentVertexIterator {
    HalfEdge he;
    bool justStarted = true;

    VertexAdjacentVertexIterator& operator++() { he = he.nextIncomingHalfEdge(); justStarted = false; return *this; }

    bool operator !=(VertexAdjacentVertexIterator const& other) const { return justStarted || he.idx != other.he.idx; }

    Vertex operator*() const { return he.tail(); }
};

template<class T>
struct ElementIterator {
    T e;

    ElementIterator& operator++() { ++e.idx; return *this; }

    ElementIterator& operator++() requires std::is_same_v<T, Corner> {
        do { ++e.idx; } while(e.idx < e.mesh->halfEdgeCount() && !isCorner());
        return *this;
    }

    bool isCorner() requires std::is_same_v<T, Corner> { return e.face().idx != Invalid; }

    bool operator !=(ElementIterator const& other) const { return e.idx != other.e.idx; }

    T operator*() const { return e; }
};

}

