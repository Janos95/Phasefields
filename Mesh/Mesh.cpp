//
// Created by janos on 8/17/20.
//

#include "Mesh.h"
#include "MeshElements.h"

#include <Corrade/Containers/StridedArrayView.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Utility/MurmurHash2.h>

#include <Magnum/Trade/MeshData.h>
#include <Magnum/MeshTools/GenerateNormals.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/Angle.h>

#include <Magnum/GL/Buffer.h>

#include <unordered_map>
#include <Magnum/Math/Color.h>

namespace Phasefield {

using namespace Magnum::Math::Literals;

struct PairHash {
    size_t operator()(std::pair<size_t, size_t> const& x) const {
        constexpr static size_t size = sizeof(std::pair<size_t, size_t>);
        static_assert(size == 2*sizeof(size_t));
        const char* pc = (char*)&x;
        return *reinterpret_cast<size_t const*>(Cr::Utility::MurmurHash2{}(pc, size).byteArray());
    }
};

inline Debug& operator<<(Debug& debug, const Implementation::HalfEdge& he) {
    char buffer[100];
    sprintf(buffer, "tail : %zu, next : %zu, face : %zu\n", he.tail, he.next, he.face);
    return debug << buffer;
}

void Mesh::setFromData(const Mg::Trade::MeshData& meshData) {
    CORRADE_ASSERT(meshData.isIndexed(), "Mesh needs to be indexed",);
    m_vertexCount = meshData.vertexCount();
    m_cornerCount = meshData.indexCount();
    m_faceCount = meshData.indexCount() / 3;
    m_dualEdgeCount = 0;
    m_halfEdgeCount = 0;

    arrayResize(m_attributes, NoInit, m_vertexCount);
    arrayResize(m_indices, NoInit, m_faceCount*3);

    /* setup index array */
    Cr::Utility::copy(arrayCast<const UnsignedInt>(meshData.indexData()), m_indices);

    /* setup attribute array */
    bool hasNormals = meshData.hasAttribute(Mg::Trade::MeshAttribute::Normal);
    for(size_t i = 0; i < m_vertexCount; ++i) {
        Implementation::Attributes attributes;
        Vector3 const& p = meshData.attribute<Vector3>(Mg::Trade::MeshAttribute::Position)[i];
        attributes.position = p;
        if(hasNormals) {
            Vector3 const& n = meshData.attribute<Vector3>(Mg::Trade::MeshAttribute::Normal)[i];
            attributes.normal = n;
        }
        m_attributes[i] = attributes;
    }
    Mg::MeshTools::generateSmoothNormalsInto(m_indices, positions(), normals());


    /* setup connectivity information */
    arrayResize(m_halfEdges, 0);
    arrayReserve(m_halfEdges, m_indices.size());
    arrayResize(m_vertexHalfEdge, NoInit, m_vertexCount);
    arrayResize(m_faceHalfEdge, NoInit, m_faceCount);
    arrayResize(degree, ValueInit, m_vertexCount);
    std::unordered_map<std::pair<size_t, size_t>, size_t, PairHash> map;

    //Debug{} << m_indices;
    //Debug{} << positions();

    for(size_t faceIdx = 0; faceIdx < m_faceCount; ++faceIdx) {
        size_t heIdx[3];
        for(size_t i = 0; i < 3; ++i) {
            size_t u = m_indices[3*faceIdx + i];
            size_t v = m_indices[3*faceIdx + (i + 1)%3];

            auto it = map.find(std::pair{v, u});
            if (it != map.end()) {
                heIdx[i] = it->second ^ size_t{1};
                m_halfEdges[heIdx[i]].face = faceIdx;
                ++m_dualEdgeCount; /* incident face was already handled, so there must be a dual edge here */
            } else {
                Implementation::HalfEdge he, twin;
                he.face = faceIdx;
                he.tail = u;
                twin.tail = v;

                arrayAppend(m_halfEdges, {he, twin});
                heIdx[i] = m_halfEdges.size() - 2;
                map[std::pair{u, v}] = heIdx[i];
            }

            m_vertexHalfEdge[u] = heIdx[i]; /* doesn't matter which one we take */
            ++degree[u];
        }

        /* hook up face to first half edge */
        m_faceHalfEdge[faceIdx] = heIdx[0];

        for(size_t i = 0; i < 3; ++i)
            m_halfEdges[heIdx[i]].next = heIdx[(i + 1)%3];
    }



    m_halfEdgeCount = m_halfEdges.size();
    m_edgeCount = m_halfEdgeCount/2;

    /* connect boundary loops */
    std::unordered_map<size_t, size_t> boundaryMap;

    for(size_t i = 0; i < m_halfEdgeCount; ++i) {
        auto& he = m_halfEdges[i];
        if(he.next == Invalid) {
            size_t tail = m_halfEdges[i].tail;
            boundaryMap[tail] = i;

        }
    }

    for(size_t i = 0; i < m_halfEdgeCount; ++i) {
        auto& he = m_halfEdges[i];
        if(he.next == Invalid) {
            size_t tip = m_halfEdges[i ^ 1].tail;
            auto it = boundaryMap.find(tip);
            CORRADE_ASSERT(it != boundaryMap.end(), "Mesh is not manifold",);
            he.next = it->second;
        }
    }

    update();
    CORRADE_ASSERT(isTriangularMesh(), "Mesh is not triangular", );
    CORRADE_ASSERT(checkDualGraph(), "Dual Graph is not correct", );
}

Mesh::Mesh(Mg::Trade::MeshData const& meshData) { setFromData(meshData); }


StridedArrayView1D<const Vector3> Mesh::positions() const {
    StridedArrayView1D<const Implementation::Attributes> view{m_attributes};
    return view.slice<const Vector3>(&Implementation::Attributes::position);
}

StridedArrayView1D<Vector3> Mesh::positions() {
    return stridedArrayView(m_attributes).slice(&Implementation::Attributes::position);
}

StridedArrayView1D<const Vector3> Mesh::normals() const {
    StridedArrayView1D<const Implementation::Attributes> view{m_attributes};
    return view.slice<const Vector3>(&Implementation::Attributes::normal);
}

StridedArrayView1D<Vector3> Mesh::normals() {
    return stridedArrayView(m_attributes).slice(&Implementation::Attributes::normal);
}

StridedArrayView1D<const Vector2> Mesh::textureCoordinates() const {
    StridedArrayView1D<const Implementation::Attributes> view{m_attributes};
    return view.slice<const Vector2>(&Implementation::Attributes::uv);
}

StridedArrayView1D<Vector2> Mesh::textureCoordinates() {
    return stridedArrayView(m_attributes).slice(&Implementation::Attributes::uv);
}

namespace {
auto sliceScalarsConst(ArrayView<const Implementation::Attributes> data) {
    auto uvs = stridedArrayView(data).slice<const Vector2>(&Implementation::Attributes::uv);
    return arrayCast<2, const float>(uvs).slice<1>();
}
auto sliceScalars(ArrayView<Implementation::Attributes> data) {
    auto uvs = stridedArrayView(data).slice(&Implementation::Attributes::uv);
    return arrayCast<2, float>(uvs).slice<1>();
}
}

StridedArrayView1D<const Float> Mesh::scalars() const { return sliceScalarsConst(m_attributes); }

StridedArrayView1D<Float> Mesh::scalars() { return sliceScalars(m_attributes); }

StridedArrayView1D<const Color4> Mesh::colors() const {
    StridedArrayView1D<const Implementation::Attributes> view{m_attributes};
    return view.slice<const Color4>(&Implementation::Attributes::color);
}

StridedArrayView1D<Color4> Mesh::colors() {
    return stridedArrayView(m_attributes).slice(&Implementation::Attributes::color);
}

VertexSet Mesh::vertices() { return {{0, this}, {m_vertexCount, this}}; }

FaceSet Mesh::faces() { return {{0, this}, {m_faceCount, this}}; }

EdgeSet Mesh::edges() { return {{0, this}, {m_edgeCount, this}}; }

HalfEdgeSet Mesh::halfEdges() { return {{0, this}, {m_halfEdgeCount, this}}; }

DualEdgeSet Mesh::dualEdges() {
    DualEdgeIterator b{0, this};
    if(m_edgeCount && !b.e.isValid()) ++b;
    return {b, {m_edgeCount, this}};
}

CornerSet Mesh::corners() {
    CornerIterator b{0, this};
    if(m_halfEdgeCount && !b.e.isValid()) ++b;
    return {b, {m_halfEdgeCount, this}};
}

void Mesh::uploadVertexBuffer(Mg::GL::Buffer& buffer) const { buffer.setData(m_attributes); }

void Mesh::uploadIndexBuffer(Mg::GL::Buffer& buffer) const { buffer.setData(m_indices); }

size_t Mesh::faceCount() const { return m_faceCount; }

size_t Mesh::vertexCount() const { return m_vertexCount; }

size_t Mesh::edgeCount() const { return m_edgeCount; }

size_t Mesh::indexCount() const { return m_indices.size(); }

size_t Mesh::halfEdgeCount() const { return m_halfEdgeCount; }

size_t Mesh::angleCount() const { return m_cornerCount; }

size_t Mesh::dualEdgeCount() const { return m_dualEdgeCount; }

size_t Mesh::lookUpFeature(const char* lookUpName) {
    for(size_t i = 0; i < m_features.size(); ++i) {
        if(strcmp(m_features[i]->name(), lookUpName) == 0)
            return i;
    }
    return Invalid;
}

void Mesh::update() {
    for(MeshFeature* feature : m_features) {
        Debug{} << "Updating feature" << feature->name();
        feature->update();
    }
}

Mg::Trade::MeshData Mesh::meshDataView() const {
    Mg::Trade::MeshIndexData indexData{m_indices};
    Array<Mg::Trade::MeshAttributeData> m_attributeDescription{InPlaceInit, {
            Mg::Trade::MeshAttributeData{Mg::Trade::MeshAttribute::Position, positions()},
            Mg::Trade::MeshAttributeData{Mg::Trade::MeshAttribute::Normal, normals()},
            Mg::Trade::MeshAttributeData{Mg::Trade::MeshAttribute::TextureCoordinates, textureCoordinates()},
            Mg::Trade::MeshAttributeData{Mg::Trade::MeshAttribute::Color, colors()},
    }};

    return Mg::Trade::MeshData{Mg::MeshPrimitive::Triangles, {},
                               m_indices,
                               indexData,
                               {},
                               m_attributes,
                               std::move(m_attributeDescription)};
}

Mesh::~Mesh() {
    for(MeshFeature* feature : m_features) {
        if(feature->ownedByMesh()) delete feature;
    }
}

ArrayView<const Vector3ui> Mesh::triangels() const { return arrayCast<const Vector3ui>(m_indices); }

bool Mesh::isTriangularMesh() {
    for(Face face : faces()) {
        HalfEdge he = face.halfEdge();
        HalfEdge it = he;
        for(size_t i = 0; i < 3; ++i)
            it = it.next();
        if(he != it) return false;
    }
    return true;
}

bool Mesh::checkDualGraph() {

    auto haveCommonEdge = [](Face a, Face b) {
        for(Edge e1 : a.edges())
        for(Edge e2 : b.edges())
            if(e1 == e2) return true;
        return false;
    };

    for(Face face : faces()) {
        for(DualEdge de : face.dualEdges()) {
            Face other = de.otherFace(face);
            if(!haveCommonEdge(face, other)) return false;
        }
    }
    return true;
}

Color4& Mesh::color(Vertex const& v) { return colors()[v.idx]; }

Color4 const& Mesh::color(Vertex const& v) const { return colors()[v.idx]; }

Float& Mesh::scalar(const Vertex& v) { return scalars()[v.idx]; }

}
