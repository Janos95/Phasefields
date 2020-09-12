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
};

struct Vertex {

    size_t idx;
    Mesh* mesh;

    auto operator<=>(Vertex const& other) const { return idx <=> other.idx; }
    bool operator==(Vertex const& other) const = default;

    [[nodiscard]] Float gaussianCurvature() const;

    [[nodiscard]] Vector3& normal() const;

    [[nodiscard]] Vector3& position() const;

    [[nodiscard]] Float& scalar() const;

    [[nodiscard]] HalfEdge halfEdge() const;

    [[nodiscard]] VertexCornerRange corners() const;

    [[nodiscard]] VertexAdjacentVertexRange adjacentVertices() const;

    [[nodiscard]] IncomingHalfEdgeRange incomingHalfEdges() const;

    [[nodiscard]] OutgoingHalfEdgeRange outgoingHalfEdges() const;

    //IncidentEdgeRange incidentEdges() const;
    //IncidentFaceRange incidentFaces() const;

};

struct Face {

    size_t idx;
    Mesh* mesh;

    auto operator<=>(Face const& other) const { return idx <=> other.idx; }
    bool operator==(Face const& other) const = default;

    //FaceIncidentVertexRange vertices();
    //FaceIncidentFaceRange incidentFaces() const;

    [[nodiscard]] HalfEdge halfEdge() const;

    [[nodiscard]] FaceEdgeRange edges() const;

    [[nodiscard]] FaceCornerRange corners() const;

    [[nodiscard]] FaceAdjacentFaceRange adjacentFaces() const;

    [[nodiscard]] FaceVertexRange vertices() const;

    [[nodiscard]] FaceHalfEdgeRange halfEdges() const;

    [[nodiscard]] bool isValid() const;

    [[nodiscard]] double area() const;

    [[nodiscard]] Vector3d normal() const;
};

struct Edge {

    size_t idx;
    Mesh* mesh;

    auto operator<=>(Edge const& other) const { return idx <=> other.idx; }
    bool operator==(Edge const& other) const = default;

    //bool isManifold();
    Vertex vertex1() const;

    Vertex vertex2() const;

    HalfEdge halfEdge() const;
};


struct Corner {

    size_t idx;
    Mesh* mesh;

    auto operator<=>(Corner const& other) const { return idx <=> other.idx; }
    bool operator==(Corner const& other) const = default;

    [[nodiscard]] HalfEdge side1() const;

    [[nodiscard]] HalfEdge side2() const;

    [[nodiscard]] Vertex vertex() const;

    [[nodiscard]] Face face() const;

    [[nodiscard]] Rad angle() const;

    [[nodiscard]] HalfEdge halfEdge() const;

};

Debug& operator<<(Debug& debug, Vertex const& v);
Debug& operator<<(Debug& debug, Face const& f);
Debug& operator<<(Debug& debug, HalfEdge const& he);
Debug& operator<<(Debug& debug, Edge const& e);
Debug& operator<<(Debug& debug, Corner const& c);

template<class E>
struct FaceCirculationIterator {
    HalfEdge he;
    bool justStarted = true;

    FaceCirculationIterator& operator++() { he = he.next(); justStarted= false; return *this;}

    bool operator !=(FaceCirculationIterator const& other) const { return justStarted || he != other.he;  }

    E operator*() const requires std::is_same_v<E, Edge> { return he.edge(); }
    E operator*() const requires std::is_same_v<E, Face> { return he.twin().face(); }
    E operator*() const requires std::is_same_v<E, HalfEdge> { return he; }
    E operator*() const requires std::is_same_v<E, Vertex> { return he.tail(); }
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

    ElementIterator& operator++() { ++e.idx; return *this; }

    ElementIterator& operator++() requires std::is_same_v<T, Corner> {
        do { ++e.idx; } while(e.idx < e.mesh->halfEdgeCount() && !isCorner());
        return *this;
    }

    bool isCorner() requires std::is_same_v<T, Corner> { return e.face().isValid(); }

    bool operator !=(ElementIterator const& other) const { return e != other.e; }

    T operator*() const { return e; }
};

}

