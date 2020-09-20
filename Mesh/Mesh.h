//
// Created by janos on 8/17/20.
//

#pragma once

#include "Surface.h"
#include "Types.h"
#include "SparseMatrix.h"
#include "MeshFeature.h"
#include "FEM.h"
#include "MeshElements.h"

#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/Pointer.h>

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
 * wrapper around an Array and ArrayView which support type safe indexing
 * and also check for out of bounds array access which can be disabled by
 * defining CORRADE_NO_ASSERT
 * @tparam N
 * @tparam T
 */
template<class E, class T>
class MeshData : public Array<T> {
public:

    T& operator[](E n) {
        CORRADE_CONSTEXPR_ASSERT(n.idx < this->size(), "Index Out of Bounds");
        return (*this)[n.idx];
    }

    T const& operator[](E n) const {
        CORRADE_CONSTEXPR_ASSERT(n.idx < this->size(), "Index Out of Bounds");
        return (*this)[n.idx];
    }

    using Array<T>::Array;
};

template<class E, class T>
class MeshDataView : public ArrayView<T> {
public:

    T& operator[](E n) {
        CORRADE_CONSTEXPR_ASSERT(n.idx < this->size(), "Index Out of Bounds");
        return (*this)[n.idx];
    }

    T const& operator[](E n) const {
        CORRADE_CONSTEXPR_ASSERT(n.idx < this->size(), "Index Out of Bounds");
        return (*this)[n.idx];
    }

    using ArrayView<T>::ArrayView;
};

template<class T> using VertexData = MeshData<Vertex, T>;
template<class T> using VertexDataView = MeshDataView<Vertex, T>;

template<class T> using FaceData = MeshData<Face, T>;
template<class T> using FaceDataView = MeshDataView<Face, T>;

template<class T> using EdgeData = MeshData<Edge, T>;
template<class T> using EdgeDataView = MeshDataView<Edge, T>;

template<class T> using HalfEdgeData = MeshData<HalfEdge, T>;
template<class T> using HalfEdgeDataView = MeshDataView<HalfEdge, T>;

template<class T> using CornerData = MeshData<Corner, T>;
template<class T> using CornerDataView = MeshDataView<Corner, T>;

template<class T> using DualEdgeData = MeshData<DualEdge, T>;
template<class T> using DualEdgeDataView = MeshDataView<DualEdge, T>;

class Mesh {
public:

    friend Vertex;
    friend HalfEdge;
    //friend Edge;
    friend Face;
    friend MeshFeature;

    explicit Mesh() = default;

    explicit Mesh(Mg::Trade::MeshData const& meshData);

    void setFromData(Mg::Trade::MeshData const& meshData);

    template<class T, class... Args>
    void requireFeature(Args&&... args) {
        auto feature = pointer<T>((Args&&)args...);
        size_t idx = lookUpFeature(feature->name());
        if(idx == Invalid) {
            //Debug{} << "Requiring feature" << feature->name();
            feature->update();
            arrayAppend(m_features, InPlaceInit, std::move(feature));
        }
    }

    void requireEdgeLengths() { requireFeature<EdgeLengthFeature>(*this); }

    void requireFaceAreas() { requireFeature<FaceAreaFeature>(*this); }

    void requireAngles() { requireFeature<AngleFeature>(*this); }

    void requireGaussianCurvature() {
        requireAngles();
        requireFeature<GaussianCurvatureFeature>(*this);
    }

    void requireGradientOperator() {
        requireFaceAreas();
        requireFeature<GradientFeature>(*this);
    }

    void requireMassMatrix() {
        requireFaceAreas();
        requireFeature<MassMatrixFeature>(*this);
    }

    void requireIntegralOperator() {
        requireFaceAreas();
        requireFeature<IntegralOperatorFeature>(*this);
    }

    [[nodiscard]] size_t vertexCount() const;

    [[nodiscard]] size_t faceCount() const;

    [[nodiscard]] size_t edgeCount() const;

    [[nodiscard]] size_t indexCount() const;

    [[nodiscard]] size_t halfEdgeCount() const;

    [[nodiscard]] size_t angleCount() const;

    [[nodiscard]] size_t dualEdgeCount() const;

    VertexSet vertices();

    FaceSet faces();

    EdgeSet edges();

    HalfEdgeSet halfEdges();

    CornerSet corners();

    DualEdgeSet dualEdges();

    void uploadVertexBuffer(Mg::GL::Buffer& vertexBuffer) const;

    void uploadScalars(Mg::GL::Buffer& vertexBuffer) const;

    void uploadIndexBuffer(Mg::GL::Buffer& indexBuffer) const;

    [[nodiscard]] StridedArrayView1D<const Vector3> positions() const;

    [[nodiscard]] StridedArrayView1D<Vector3> positions();

    [[nodiscard]] StridedArrayView1D<const Vector3> normals() const;

    [[nodiscard]] StridedArrayView1D<Vector3> normals();

    [[nodiscard]] StridedArrayView1D<const Float> scalars() const;

    [[nodiscard]] StridedArrayView1D<Float> scalars();

    [[nodiscard]] StridedArrayView1D<const Color4> colors() const;

    [[nodiscard]] StridedArrayView1D<Color4> colors();

    void update();

    EdgeData<double> edgeLength;

    CornerData<Rad> angle;

    FaceData<double> faceArea;
    FaceData<Vector3d> faceNormal;
    FaceData<double> faceDiameter;

    VertexData<double> gaussianCurvature;
    VertexData<size_t> degree;

    /* Fem operators */
    VertexData<double> integral;
    Array<Triplet> massMatrix;
    //Array<Triplet> stiffnessMatrix;
    HalfEdgeData<Vector3d> gradient;

    double surfaceArea;

protected:

    size_t lookUpFeature(const char*);

    bool m_requireEdgeLength;
    bool m_requireAngle;
    bool m_requiresFaceArea;
    bool m_requiresGaussianCurvature;
    bool m_requireDegree;
    bool m_requireIntegralOperator;

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
    size_t m_dualEdgeCount = 0;

    Array<Pointer<MeshFeature>> m_features;
};

}
