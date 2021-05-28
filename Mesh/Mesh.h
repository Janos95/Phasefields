//
// Created by janos on 8/17/20.
//

#pragma once

#include "Surface.h"
#include "Types.h"
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
 * and also check for incorrect array access which can be disabled by
 * defining CORRADE_NO_ASSERT
 * @tparam N
 * @tparam T
 */


template<class E, class T>
class MeshDataView : public ArrayView<T> {
public:

    /* implicit */ MeshDataView(ArrayView<T> v) : ArrayView<T>(v) {}

    T& operator[](E n) {
        CORRADE_CONSTEXPR_ASSERT(n.idx < this->size(), "Index Out of Bounds");
        return (*this)[n.idx];
    }

    T const& operator[](E n) const {
        CORRADE_CONSTEXPR_ASSERT(n.idx < this->size(), "Index Out of Bounds");
        return (*this)[n.idx];
    }

    using ArrayView<T>::ArrayView;
    operator ArrayView<T>() {
        auto* p = reinterpret_cast<ArrayView<T>*>(this);
        return *p;
    }
};

template<class E, class T>
class MeshData : public Array<T> {
public:

    T& operator[](E n) {
        CORRADE_CONSTEXPR_ASSERT(n.idx < this->size(), "Index Out of Bounds");
        if constexpr (std::is_same_v<E, Corner> || std::is_same_v<E, DualEdge>)
            CORRADE_CONSTEXPR_ASSERT(n.isValid(), "Element does not represent a valid index");
        return (*this)[n.idx];
    }

    T const& operator[](E n) const {
        CORRADE_CONSTEXPR_ASSERT(n.idx < this->size() && n, "Index Out of Bounds");
        return (*this)[n.idx];
    }

    using Array<T>::Array;

    operator MeshDataView<E, T>() {
        auto* p = reinterpret_cast<Array<T>*>(this);
        ArrayView<T> view = *p;
        return view;
    }

    static constexpr bool HasConst = std::is_const_v<T>;

    operator ArrayView<std::conditional_t<HasConst, const T, T>>() {
        auto* p = reinterpret_cast<Array<T>*>(this);
        return *p;
    }
};

class Mesh {
public:

    friend Vertex;
    friend HalfEdge;
    //friend Edge;
    friend Face;
    friend MeshFeature;

    explicit Mesh() = default;
     ~Mesh();

    explicit Mesh(Mg::Trade::MeshData const& meshData);

    /**
     * Copying is not allowed. Use the clone method instead to be more
     * explicit.
     */
    Mesh(Mesh const&) = delete;
    Mesh& operator=(Mesh const&) = delete;

    Mesh(Mesh&&) noexcept = default;
    Mesh& operator=(Mesh&&) noexcept = default;

    void setFromData(Mg::Trade::MeshData const& meshData);

    [[nodiscard]] Mg::Trade::MeshData meshDataView() const;

    template<class T>
    void checkAndAddFeature() {
        size_t idx = lookUpFeature(T::Name);
        if(idx == Invalid) {
            MeshFeature* f = new T{*this};
            f->update();
        }
    }

    void requireEdgeLengths() { checkAndAddFeature<EdgeLengthFeature>(); }
    void requireFaceInformation() { checkAndAddFeature<FaceInformationFeature>(); }
    void requireAngles() { checkAndAddFeature<AngleFeature>(); }
    void requireBoundaryInformation() { checkAndAddFeature<BoundaryInformation>(); }
    void requireGaussianCurvature() {
        requireAngles();
        requireIntegralOperator();
        requireBoundaryInformation();
        checkAndAddFeature<GaussianCurvatureFeature>();
    }
    void requireGradientOperator() {
        requireFaceInformation();
        checkAndAddFeature<GradientFeature>();
    }
    void requireMassMatrix() {
        requireFaceInformation();
        checkAndAddFeature<MassMatrixFeature>();
    }
    void requireStiffnessMatrix() {
        requireBoundaryInformation();
        requireFaceInformation();
        requireIntegralOperator();
        checkAndAddFeature<StiffnessMatrixFeature>();
    }
    void requireIntegralOperator() {
        requireFaceInformation();
        checkAndAddFeature<IntegralOperatorFeature>();
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

    /**
     * Flips the edge, if the edge is not on the boundary loop.
     * The Edge is NOT invalidated, so it can be safely reused.
     */
    void flip(Edge);

    /**
     * Splits an edge. This does not invalidate the edge, but
     * the edge now references one of edges resulting from the split.
     */
    HalfEdge split(Edge);

    HalfEdge insertVertexAlongEdge(Edge);

    void connectVertices(HalfEdge, HalfEdge);

    void collapse(Edge);

    void loopSubdivide();

    void catmullClark();

    void isotropicRemeshing();

    void triangulate();

    void triangulate(Face f);

    void generateFaceList();

    void computeVertexNormals();

    //Mesh clone();

    Vertex makeVertex();

    Face makeFace();

    void compress();

    void uploadVertexBuffer(Mg::GL::Buffer& vertexBuffer) const;

    void uploadIndexBuffer(Mg::GL::Buffer& indexBuffer) const;

    [[nodiscard]] StridedArrayView1D<const Vector3> positions() const;

    [[nodiscard]] StridedArrayView1D<Vector3> positions();

    [[nodiscard]] StridedArrayView1D<const Vector3> normals() const;

    [[nodiscard]] StridedArrayView1D<Vector3> normals();

    [[nodiscard]] StridedArrayView1D<const Float> scalars() const;

    [[nodiscard]] StridedArrayView1D<Float> scalars();

    [[nodiscard]] StridedArrayView1D<const Color4> colors() const;

    [[nodiscard]] StridedArrayView1D<Color4> colors();

    [[nodiscard]] Color4& color(Vertex const&);

    [[nodiscard]] Float& scalar(Vertex const&);

    [[nodiscard]] Vector3& position(Vertex const&);

    [[nodiscard]] Vector3& normal(Vertex const&);

    [[nodiscard]] StridedArrayView1D<const Vector2> textureCoordinates() const;

    [[nodiscard]] StridedArrayView1D<Vector2> textureCoordinates();

    [[nodiscard]] ArrayView<const Vector3ui> triangels() const;

    [[nodiscard]] ArrayView<const UnsignedInt> indices() const;

    void update();

    EdgeData<double> edgeLength;

    CornerData<Radd> angle;

    FaceData<double> faceArea;
    FaceData<Vector3d> faceNormal;
    FaceData<double> faceDiameter;
    double surfaceArea;

    VertexData<double> gaussianCurvature;
    FaceData<double> faceCurvature;
    VertexData<size_t> degree;

    VertexData<bool> isOnBoundary;
    HalfEdgeData<double> stiffnessElements;

    /* Fem operators */
    VertexData<double> integral;
    HalfEdgeData<Vector3d> gradient;

    bool isTriangularMesh();
    bool checkDualGraph();

protected:

    size_t lookUpFeature(const char*);

    Array<Implementation::Attributes> m_attributes;
    Array<UnsignedInt> m_indices;

    HalfEdgeData<Implementation::HalfEdge> m_halfEdges;
    FaceData<UnsignedInt> m_faceHalfEdge;
    VertexData<UnsignedInt> m_vertexHalfEdge;

    size_t m_halfEdgeCount = 0;
    size_t m_faceCount = 0;
    size_t m_vertexCount = 0;

    size_t m_dualEdgeCount = 0;

    Array<MeshFeature*> m_features;
};

}
