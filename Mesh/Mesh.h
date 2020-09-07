//
// Created by janos on 8/17/20.
//

#pragma once

#include "Surface.h"
#include "Types.h"

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/LinkedList.h>

#include <Magnum/Trade/Trade.h>
#include <Magnum/GL/GL.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Math/Color.h>

namespace Phasefield {

namespace Mg = Magnum;
namespace Cr = Corrade;

namespace Implementation {

struct HalfEdge {
    size_t next = Invalid;
    size_t face = Invalid;
    size_t tail = Invalid;
};

struct Attributes {
    Vector3 position;
    Vector3 normal;
    Vector2 uv;
    Color4 color;
};

struct HalfEdgeIncidency {
    size_t in;
    size_t out;
};

}

/**
 * wrapper around an Array<T> which can be indexed with an associated Navigator N
 * @tparam N
 * @tparam T
 */
template<class N, class T>
class MeshData : public Array<T> {
public:
    T& operator[](N n) { return (*this)[n.idx]; }
    T const& operator[](N n) const { return (*this)[n.idx]; }

    using Array<T>::Array;
};

template<class T>
using VertexData = MeshData<Vertex, T>;

template<class T>
using FaceData = MeshData<Face, T>;

template<class T>
using EdgeData = MeshData<Edge, T>;

template<class T>
using HalfEdgeData = MeshData<HalfEdge, T>;

template<class T>
using CornerData = MeshData<Corner, T>;


class Mesh {
public:

    //class Feature : public LinkedListItem<Feature> {
    //public:
    //    Feature(Mesh& mesh) : m_mesh(&mesh) {}
    //    virtual ~Feature() { m_mesh->features.erase()}
    //    virtual void update() = 0;
    //private:
    //    Mesh* m_mesh = nullptr;

    //};

    friend Vertex;
    friend HalfEdge;
    //friend Corner;
    //friend Edge;
    friend Face;
    //friend Feature;

    explicit Mesh() = default;

    explicit Mesh(Mg::Trade::MeshData const& meshData);

    void setFromData(Mg::Trade::MeshData const& meshData);

    void requireEdgeLengths();

    void requireFaceAreas();

    void requireAngles();

    void requireGaussianCurvature();

    [[nodiscard]] size_t vertexCount() const;

    [[nodiscard]] size_t faceCount() const;

    [[nodiscard]] size_t edgeCount() const;

    [[nodiscard]] size_t indexCount() const;

    [[nodiscard]] size_t halfEdgeCount() const;

    VertexSet vertices();

    FaceSet faces();

    EdgeSet edges();

    HalfEdgeSet halfEdges();

    CornerSet corners();

    void uploadVertexBuffer(Mg::GL::Buffer& vertexBuffer) const;

    void uploadIndexBuffer(Mg::GL::Buffer& indexBuffer) const;

    StridedArrayView1D<Vector3> positions();

    StridedArrayView1D<Vector3> normals();

    StridedArrayView1D<Float> scalars();

    StridedArrayView1D<Color4> colors();

    EdgeData<float> edgeLength;
    VertexData<float> gaussianCurvature;
    CornerData<Rad> angle;
    VertexData<size_t> degree;
    FaceData<float> faceArea;

private:

    Array<Implementation::Attributes> m_attributes;
    Array<UnsignedInt> m_indices;

    Array<Implementation::HalfEdge> m_halfEdges;
    Array<UnsignedInt> m_faceHalfEdge;
    Array<UnsignedInt> m_vertexHalfEdge;

    Array<size_t> m_outgoingHalfEdges;

    size_t m_faceCount = 0;
    size_t m_vertexCount = 0;
    size_t m_edgeCount = 0;
    size_t m_halfEdgeCount = 0;
    size_t m_cornerCount = 0;

    //LinkedList<Feature> features;
};

}
