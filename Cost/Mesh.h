//
// Created by janos on 8/17/20.
//

#ifndef PHASEFIELD_MESH_H
#define PHASEFIELD_MESH_H

#include "Types.h"

#include <Corrade/Containers/Array.h>
#include <Magnum/Trade/Trade.h>

namespace Phasefield {

namespace Mg = Magnum;
namespace Cr = Corrade;

namespace Implementation {
struct HalfEdge {
    static constexpr UnsignedInt Invalid = ~0u;

    UnsignedInt next;
    UnsignedInt face;
    UnsignedInt opposite;
    UnsignedInt vertex;
};

struct Attributes {
    Vector3 position;
    Vector3 normal;
    Float scalar;
};

struct HalfEdgeIncidency {
    UnsignedInt nextOut;
    UnsignedInt nextIn;
};

}

class Mesh {
public:

    struct Vertex;
    struct Edge;
    struct HalfEdge;
    struct Face;

    template<class T>
    struct IncidentHalfEdgeIterator {
        Int idx;
        Array<UnsignedInt>* next;

        IncidentHalfEdgeIterator& operator++() {
            idx = (*next)[idx];
            return *this;
        }

        T operator*();

    };

    template<class T>
    struct IncidentHalfEdgeRange {

        IncidentHalfEdgeIterator<T> begin() {

        }

        IncidentHalfEdgeIterator<T> end() {

        }
    };

    template<class T>
    struct ElementSetIterator {
        UnsignedInt idx;
        Mesh* mesh;

        T operator*() { return T{this, idx}; }
        auto operator<=>(const ElementSetIterator& other) const { return idx <=> other.idx; }
    };

    template<class T>
    struct ElementSet {
        using Iterator = ElementSetIterator<T>;
        UnsignedInt elementCount;
        Mesh *mesh;

        Iterator begin() { return {0, mesh}; }
        Iterator end() { return {elementCount, mesh}; }
    };


    using IncidentEdgeRange = IncidentHalfEdgeRange<Edge>;
    using AdjacentVertexRange = IncidentHalfEdgeRange<Vertex>;
    using IncidentFaceRange = IncidentHalfEdgeRange<Face>;
    using IncomingHalfEdgeRange = IncidentHalfEdgeRange<HalfEdge>;
    using OutgoingHalfEdgeRange = IncidentHalfEdgeRange<HalfEdge>;

    explicit Mesh(Mg::Trade::MeshData meshData);

    struct HalfEdge {

        HalfEdge next() const;
        Vertex tip() const;
        Vertex tail() const;
        Face face() const;

        Mg::UnsignedInt idx;
        Mesh* mesh;
    };

    struct Vertex {

        Float& gaussianCurvature();
        Vector3& normal();
        Vector3& position();

        HalfEdge halfEdge();

        AdjacentVertexRange adjacentVertices() const;
        IncomingHalfEdgeRange incomingHalfedges() const;
        OutgoingHalfEdgeRange outgoingHalfedges() const;
        IncidentEdgeRange incidentEdges() const;
        IncidentFaceRange incidentFaces() const;

        Mg::UnsignedInt idx;
        Mesh* mesh;
    };

    struct Face {

        FaceIncidentVertexRange vertices();
        FaceIncidentFaceRange incidentFaces() const;

        Vector3& a();
        Vector3& b();
        Vector3& c();

        Mg::UnsignedInt idx;
        Mesh* mesh;
    };

    struct Edge {

        bool isManifold();

        Vertex tip();
        Vertex tail();
        HalfEdge halfEdge();

        UnsignedInt idx;
        Mesh* mesh;
    };

    struct AdjacentVertexIterator {
        HalfEdge he;

        auto operator<=>(const AdjacentVertexIterator&) const = default;
    };

    using FaceSet = ElementSet<Face>;
    using VertexSet = ElementSet<Vertex>;
    using EdgeSet = ElementSet<Edge>;
    using HalfEdgeSet = ElementSet<HalfEdge>;

    bool isManifold();

    void requireEdgeLengths();

    void requireInternalAngles();

    void requireGaussianCurvature();

    UnsignedInt vertexCount() const;

    UnsignedInt faceCount() const;

    VertexSet vertices() { return VertexSet{m_vertexCount, this}; }

    FaceSet faces() { return FaceSet{m_faceCount, this}; }

    EdgeSet edges() { return EdgeSet{m_edgeCount, this}; }

    HalfEdgeSet halfEdges() { return HalfEdgeSet{m_halfEdgeCount, this}; }

    StridedArrayView1D<Vector3> positions();

    StridedArrayView1D<Vector3> normals();

    StridedArrayView1D<Float> scalarField();

    Array<Float>& gaussianCurvature();

    Array<Float>& angles();

private:

    Array<Float> m_edgeLengths;
    Array<Float> m_gaussianCurvature;

    Array<UnsignedInt> m_indices;
    Array<Rad> m_angles;

    Array<Implementation::Attributes> m_attributes;

    Array<Implementation::HalfEdge> m_halfEdges;
    Array<UnsignedInt> m_faceHalfEdge;
    Array<UnsignedInt> m_edgeHalfEdge;
    Array<UnsignedInt> m_vertexHalfEdge;

    /* this are not strictly necessary */
    Array<Implementation::HalfEdgeIncidency> m_outgoingHalfEdges;

    UnsignedInt m_faceCount;
    UnsignedInt m_vertexCount;
    UnsignedInt m_edgeCount;
    UnsignedInt m_halfEdgeCount;
};


}


#endif //PHASEFIELD_MESH_H
