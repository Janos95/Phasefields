//
// Created by janos on 8/17/20.
//

#include "Mesh.h"
#include "MeshElements.h"
#include "Algorithms.h"

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
    size_t vertexCount = meshData.vertexCount();
    m_dualEdgeCount = 0;
    size_t faceCount = meshData.indexCount()/3;

    arrayResize(m_attributes, NoInit, vertexCount);
    arrayResize(m_indices, NoInit, meshData.indexCount());

    /* setup index array */
    Cr::Utility::copy(arrayCast<const UnsignedInt>(meshData.indexData()), m_indices);

    /* setup attribute array */
    bool hasNormals = meshData.hasAttribute(Mg::Trade::MeshAttribute::Normal);
    for(size_t i = 0; i < vertexCount; ++i) {
        Implementation::Attributes attributes;
        Vector3 const& p = meshData.attribute<Vector3>(Mg::Trade::MeshAttribute::Position)[i];
        attributes.position = p;
        if(hasNormals) {
            Vector3 const& n = meshData.attribute<Vector3>(Mg::Trade::MeshAttribute::Normal)[i];
            attributes.normal = n;
        }
        m_attributes[i] = attributes;
    }


    /* setup connectivity information */
    arrayResize(m_halfEdges, 0);
    arrayReserve(m_halfEdges, m_indices.size());
    arrayResize(m_vertexHalfEdge, NoInit, vertexCount);
    arrayResize(m_faceHalfEdge, NoInit, faceCount);
    arrayResize(degree, ValueInit, vertexCount);
    std::unordered_map<std::pair<size_t, size_t>, size_t, PairHash> map;

    //Debug{} << m_indices;
    //Debug{} << positions();

    for(size_t faceIdx = 0; faceIdx < faceCount; ++faceIdx) {
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

    /* connect boundary loops */
    std::unordered_map<size_t, size_t> boundaryMap;

    for(size_t i = 0; i < halfEdgeCount(); ++i) {
        auto& he = m_halfEdges[i];
        if(he.next == Invalid) {
            size_t tail = m_halfEdges[i].tail;
            boundaryMap[tail] = i;

        }
    }

    for(size_t i = 0; i < halfEdgeCount(); ++i) {
        auto& he = m_halfEdges[i];
        if(he.next == Invalid) {
            size_t tip = m_halfEdges[i ^ 1].tail;
            auto it = boundaryMap.find(tip);
            CORRADE_ASSERT(it != boundaryMap.end(), "Mesh is not manifold",);
            he.next = it->second;
        }
    }

    if(!hasNormals)
        computeVertexNormals();

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

VertexSet Mesh::vertices() { return {{0, this}, {vertexCount(), this}}; }

FaceSet Mesh::faces() { return {{0, this}, {faceCount(), this}}; }

EdgeSet Mesh::edges() { return {{0, this}, {edgeCount(), this}}; }

HalfEdgeSet Mesh::halfEdges() { return {{0, this}, {halfEdgeCount(), this}}; }

DualEdgeSet Mesh::dualEdges() {
    DualEdgeIterator b{0, this};
    if(edgeCount() && !b.e.isValid()) ++b;
    return {b, {edgeCount(), this}};
}

CornerSet Mesh::corners() {
    CornerIterator b{0, this};
    if(halfEdgeCount() && !b.e.isValid()) ++b;
    return {b, {halfEdgeCount(), this}};
}

void Mesh::uploadVertexBuffer(Mg::GL::Buffer& buffer) const { buffer.setData(m_attributes); }

void Mesh::uploadIndexBuffer(Mg::GL::Buffer& buffer) const { buffer.setData(m_indices); }

size_t Mesh::faceCount() const { return m_faceHalfEdge.size(); }

size_t Mesh::vertexCount() const { return m_attributes.size(); }

size_t Mesh::edgeCount() const { return m_halfEdges.size()/2; }

size_t Mesh::indexCount() const { return m_indices.size(); }

size_t Mesh::halfEdgeCount() const { return m_halfEdges.size(); }

size_t Mesh::angleCount() const { return m_indices.size(); }

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

Float& Mesh::scalar(const Vertex& v) { return scalars()[v.idx]; }

void Mesh::flip(Edge e) {
    if(e.onBoundaryLoop()) return;

    HalfEdge h0 = e.halfEdge();
    HalfEdge h1 = h0.next();
    HalfEdge h2 = h1.next();

    HalfEdge h3 = h0.twin();
    HalfEdge h4 = h3.next();
    HalfEdge h5 = h4.next();

    auto& h0Impl = m_halfEdges[h0.idx];
    auto& h1Impl = m_halfEdges[h1.idx];
    auto& h3Impl = m_halfEdges[h3.idx];
    auto& h4Impl = m_halfEdges[h4.idx];

    h0Impl.next = h2.idx;
    h0Impl.tail = h5.tail().idx;
    h1Impl.next = h3.idx;
    h1Impl.face = h3.face().idx;
    h3Impl.next = h5.idx;
    h3Impl.tail = h2.tail().idx;
    h4Impl.next = h0.idx;
    h4Impl.face = h0.face().idx;

    size_t face0 = h3.face().idx;
    size_t face1 = h0.face().idx;

    m_indices[face0] = h3.tail().idx;
    m_indices[face0 + 1] = h5.tail().idx;
    m_indices[face0 + 2] = h1.tail().idx;

    m_indices[face1] = h0.tail().idx;
    m_indices[face1 + 1] = h2.tail().idx;
    m_indices[face1 + 2] = h4.tail().idx;

    m_faceHalfEdge[face0] = h3.idx;
    m_faceHalfEdge[face1] = h0.idx;
}

HalfEdge Mesh::insertVertexAlongEdge(Edge e) {
    Vertex v = makeVertex();


    HalfEdge he0 = e.halfEdge();
    if(he0.onBoundaryLoop()) he0 = he0.twin();

    HalfEdge he2 = he0.next().next();
    while(he2.next() != he0) he2 = he2.next();

    /* may lie on the boundary */
    HalfEdge he3 = he0.twin();
    HalfEdge he4 = he3.next();

    Implementation::HalfEdge he6Impl;
    Implementation::HalfEdge he7Impl;
    auto& he0Impl = m_halfEdges[he0];
    auto& he3Impl = m_halfEdges[he3];
    auto& he2Impl = m_halfEdges[he2];

    size_t he6Idx = halfEdgeCount();
    size_t he7Idx = halfEdgeCount() + 1;

    /* setup new half-edges */
    he6Impl.face = he0.face().idx;
    he6Impl.next = he0.idx;
    he6Impl.tail = he0.tail().idx;

    he7Impl.face = he3.face().idx;
    he7Impl.next = he4.idx;
    he7Impl.tail = v.idx;

    /* reset old ones */
    he0Impl.tail = v.idx;
    he3Impl.next = he7Idx;
    he2Impl.next = he6Idx;

    /* set vertex to halfedge map for new vertex */
    m_vertexHalfEdge[v] = he0.idx;

    arrayAppend(m_halfEdges, {he6Impl, he7Impl});

    return he0;
}

void Mesh::connectVertices(HalfEdge heA, HalfEdge heB) {
    //Face fA = heA.face();
    //Face fB{faceCount(), this};

    //// Faces
    //m_faceHalfEdge[fA] = heANew.getIndex();
    //m_faceHalfEdge[fB] = heBNew.getIndex();

    //// Halfedges
    //heNextArr[heANew.getIndex()] = heB.getIndex();
    //heVertexArr[heANew.getIndex()] = vA.getIndex();
    //heFaceArr[heANew.getIndex()] = fA.getIndex();

    //heNextArr[heBNew.getIndex()] = heA.getIndex();
    //heVertexArr[heBNew.getIndex()] = vB.getIndex();
    //heFaceArr[heBNew.getIndex()] = fB.getIndex();

    //heNextArr[heAPrev.getIndex()] = heANew.getIndex();
    //heNextArr[heBPrev.getIndex()] = heBNew.getIndex();

    //// Set all other new .face pointers to fB
    //HalfEdge he = heA;
    //while(he != heBNew) {
    //    m_halfEdges[he].face = fB.idx;
    //    he = he.next();
    //}
}

HalfEdge Mesh::split(Edge e) {

    HalfEdge he = insertVertexAlongEdge(e);

    { // primary face
        HalfEdge heOther = he.next().next();
        connectVertices(he, heOther);
    }

    if (he.twin().isInterior()) { // secondary face
        HalfEdge heFirst = he.twin().next();
        HalfEdge heOther = heFirst.next().next();
        connectVertices(heFirst, heOther);
    }

    return he;
}

void Mesh::collapse(Edge e) {
    HalfEdge he0 = e.halfEdge();
    if(he0.onBoundaryLoop())
        he0 = he0.twin();

    Face face = he0.face();
    if(face.isTriangle()) {

    }
}

/**
 * 1. 4-1 subdivision of the mesh
 *  a) Split every edge of the mesh in any order whatsoever.
 *  b) Flip any new edge that touches a new vertex and an old vertex.
 * 2. Update Vertex positions
 */
void Mesh::loopSubdivide() {

    //auto isNewVertex = [vc = vertexCount()](Vertex v) -> bool {
    //    if(v.idx < vc)
    //        return false;
    //    return true;
    //};

    ///* an edge should be flipped if its endpoints touch a newly
    // * inserted vertex and an old vertex */
    //auto shouldFlip = [vc = vertexCount()](Edge e) -> bool {
    //    auto [min, max] = Math::minmax(e.vertex1().idx, e.vertex2().idx);
    //    if(min < vc && max >= vc) return true;
    //    else return false;
    //};

    ///* each edge will produce a new vertex */
    //VertexData<Vector3> newPositions(vertexCount() + edgeCount());

    //for(Edge e : edges()) {
    //    HalfEdge he = split(e);

    //    Vertex v1 = e.vertex1(), v2 = e.vertex2();
    //    position(v) = 0.5*(v1.position() + v2.position());
    //}

    //for(Edge e : edges()) {
    //    if(shouldFlip(e)) {
    //        flip(e);
    //    }
    //}

    //for(Vertex v : vertices()) {
    //    position(v) = newPositions[v];
    //}
}

//Mesh Mesh::clone() {
//
//    Array<Implementation::Attributes> attributes(NoInit, m_attributes.size());
//    Array<UnsignedInt> indices(NoInit, m_indices.size());
//    Array<Implementation::HalfEdge> halfEdges(NoInit, m_halfEdges.size());
//    FaceData<UnsignedInt> faceHalfEdge(NoInit, m_faceHalfEdge.size());
//    VertexData<UnsignedInt> vertexHalfEdge(NoInit, m_vertexHalfEdge.size());
//
//    Cr::Utility::copy(m_attributes, attributes);
//    Cr::Utility::copy(m_indices, indices);
//    Cr::Utility::copy(m_halfEdges, halfEdges);
//    Cr::Utility::copy(m_faceHalfEdge, faceHalfEdge);
//    Cr::Utility::copy(m_vertexHalfEdge, vertexHalfEdge);
//
//    Mesh mesh;
//    mesh.m_attributes = std::move(attributes);
//    mesh.m_indices = std::move(indices);
//    mesh.m_halfEdges = std::move(halfEdges);
//    mesh.m_faceHalfEdge = std::move(faceHalfEdge);
//    mesh.m_vertexHalfEdge = std::move(vertexHalfEdge);
//    mesh.m_dualEdgeCount = m_dualEdgeCount;
//
//    return mesh;
//}

Vertex Mesh::makeVertex() {
    arrayAppend(m_attributes, {});
    arrayAppend(m_vertexHalfEdge, ~0u);
    return {m_attributes.size(), this};
}

Face Mesh::makeFace() {
    arrayAppend(m_faceHalfEdge, ~0u);
    return Face{faceCount() - 1, this};
}

void Mesh::triangulate() {
    for(Face face : faces())
        triangulate(face);
}

void Mesh::triangulate(Face f) {
    HalfEdge he = f.halfEdge();
    HalfEdge it1 = he.next().next();
    HalfEdge it2 = it1.next();
    while(it2 != he) {
        connectVertices(it1, he);
        it1 = it2;
        it2 = it2.next();
    }
}

void Mesh::generateFaceList() {
    Array<UnsignedInt> indices;
    arrayReserve(indices, faceCount()*3);
    for(Face face : faces()) {
        for(Vertex v : face.vertices()) {
            arrayAppend(indices, UnsignedInt(v.idx));
        }
    }
    m_indices = std::move(indices);
}

void Mesh::computeVertexNormals() {
    for(Vertex v : vertices())
        normal(v) = Vector3{Math::ZeroInit};

    for(Face face : faces()) {
        Vector3 N = face.computeNormal();
        for(Corner corner : face.corners()) {
            normal(corner.vertex()) += N*float(corner.computeAngle());
        }
    }

    for(Vertex v : vertices())
        normal(v) = normal(v).normalized();
}

void Mesh::catmullClark() {

}

void Mesh::isotropicRemeshing() {

    //constexpr double upper = 4./3.;
    //constexpr double lower = 4./5.;
    //constexpr double weight = 0.2;

    ///* compute mean edge length */
    //double meanEdgeLength = 0;
    //for(Edge e : edges()) {
    //    meanEdgeLength += e.length();
    //}
    //meanEdgeLength /= edgeCount();

    ///* split long edges */
    //double thresholdToLong = upper*meanEdgeLength;
    //for(Edge e : edges()) {
    //    if(e.length() > thresholdToLong)  {
    //        split(e);
    //    }
    //}

    ///* split short edges */
    //double thresholdToShort = lower*meanEdgeLength;
    //for(Edge e : edges()) {
    //    if(e.length() < thresholdToShort)  {
    //        collapse(e);
    //    }
    //}

    ///* flip edges if it improves tatol valence deviation from optimum (6) */
    //for(Edge e : edges()) {
    //    if(e.onBoundaryLoop()) continue;

    //    size_t a1 = e.vertex1().computeDegree();
    //    size_t a2 = e.vertex2().computeDegree();

    //    HalfEdge he = e.halfEdge();
    //    size_t b1 = he.next().tip().computeDegree();
    //    size_t b2 = he.twin().next().tip().computeDegree();

    //    size_t deviation = Math::abs(a1 - 6) + Math::abs(a2 - 6) + Math::abs(b1 - 6) + Math::abs(b2 - 6);
    //    size_t deviationAfterFlip = Math::abs(a1 - 1 - 6) + Math::abs(a2 - 1 - 6) + Math::abs(b1 + 1 - 6) + Math::abs(b2 + 1 - 6);

    //    if(deviationAfterFlip < deviation)
    //        flip(e);
    //}

    ///* tangential smoothing */
    //VertexData<Vector3> smoothingDirection(NoInit, vertexCount());
    //computeVertexNormals();

    ///* compute smoothing direction */
    //for(Vertex v : vertices()) {
    //    Vector3 center{};
    //    size_t degree = 0;
    //    for(Vertex w : v.adjacentVertices()) {
    //        center += w.position();
    //        ++degree;
    //    }
    //    center /= float(degree);
    //    Vector3 dir = center - v.position();
    //    Vector3 normal = ;
    //    smoothingDirection[v] = dir - Math::dot(normal,dir)*normal;
    //}

    ///* relocate vertices */
    //for(Vertex v : vertices()) {
    //    position(v) = v.position() + weight*smoothingDirection[v];
    //}
}

void Mesh::compress() {

    size_t idx = 0;
    Array<size_t> map(NoInit, m_halfEdges.size());
    for(size_t i = 0; i < m_halfEdges.size(); ++i) {
        if(m_halfEdges[i].next != Invalid)
            map[i] = idx++;
        else
            map[i] = Invalid;
    }

    //auto removeInvalidElements = [](auto& arr, size_t& size) {
    //    constexpr static bool isHe = requires { arr[0].next == Invalid; };

    //    auto pred = [](auto& el) {
    //        if constexpr (isHe) return el.next == Invalid;
    //        else return el == Invalid;
    //    };

    //    size_t idx = 0;
    //    Array<size_t> map(NoInit, arr.size());
    //    for(size_t i = 0; i < arr.size(); ++i) {
    //        if(pred(arr[i]) != Invalid)
    //            map[i] = idx++;
    //        else
    //            map[i] = Invalid;
    //    }

    //    //for(size_t i = 0; i < arr.size(); ++i) {
    //    //    if(!pred(arr[i])) {
    //    //        if(isHe) {
    //    //            arr[]
    //    //        }
    //    //    }
    //    //}

    //    auto it = removeIf(arr, [](auto& el) { return el == Invalid; });
    //    size = arr.end() - it;
    //    arrayResize(arr, size);

    //    for(auto& x : arr) {

    //    }
    //};

    //removeInvalidElements(m_halfEdges, m_halfEdgeCount);
    //removeInvalidElements(m_vertexHalfEdge, m_vertexCount);
    //removeInvalidElements(m_faceHalfEdge, m_faceCount);
}

Vector3& Mesh::normal(const Vertex& v) { return normals()[v.idx]; }

Vector3& Mesh::position(const Vertex& v) { return positions()[v.idx]; }

ArrayView<const UnsignedInt> Mesh::indices() const { return m_indices; }

}
