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
    m_vertexCount = meshData.vertexCount();
    m_cornerCount = meshData.indexCount();
    m_faceCount = meshData.indexCount() / 3;
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
        /* hook up face to first half edge */
        m_faceHalfEdge[faceIdx] = 3*faceIdx;

        size_t heIdx[3];
        for(size_t i = 0; i < 3; ++i) {
            size_t u = m_indices[3*faceIdx + i];
            size_t v = m_indices[3*faceIdx + (i + 1)%3];

            auto it = map.find(std::pair{v, u});
            if (it != map.end()) {
                heIdx[i] = it->second ^ size_t{1};
                m_halfEdges[heIdx[i]].face = faceIdx;
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

    //for(Vertex v : vertices()) {
    //    HalfEdge he = v.halfEdge();
    //    printf("vertex %zu with half edge (%zu, %zu)\n", v.idx, he.tail().idx, he.tip().idx);
    //}


    //for(Vertex v : vertices()) {
    //    printf("Half-Edges outgoing from %zu : ", v.idx);
    //    for(HalfEdge he : v.outgoingHalfedges()) {
    //        printf(" (%zu, %zu) ", he.tail().idx, he.tip().idx);
    //    }
    //    printf("\n");
    //}

    //for(Vertex v : vertices()) {
    //    printf("Half-Edges incoming from %zu : ", v.idx);
    //    for(HalfEdge he : v.incomingHalfedges()) {
    //        printf(" (%zu, %zu) ", he.tail().idx, he.tip().idx);
    //    }
    //    printf("\n");
    //}

    //for(HalfEdge he : halfEdges()) {
    //    printf("half edge (%zu, %zu), face %zu\n", he.tail().idx, he.tip().idx, he.face().idx);
    //}

    //arrayResize(m_edgeHalfEdge, NoInit, m_edgeCount); /* shrink to correct size */
    //m_prefixSums[0] = 0;
    //for(size_t i = 1; i < m_prefixSums.size(); ++i)
    //    m_prefixSums[i] += degree[i - 1];

    //Array<size_t> currentPosition(NoInit, m_vertexCount);
    //Cr::Utility::copy({m_prefixSums.data(), m_prefixSums.size() - 1}, currentPosition);

    //for(size_t faceIdx = 0; faceIdx < m_faceCount; ++faceIdx) {
    //    for(size_t i = 0; i < 3; ++i) {
    //        size_t idx = 3*faceIdx + i;
    //        size_t next = 3*faceIdx + (i + 1)%3;
    //        m_incidentHalfEdges[currentPosition[m_indices[idx]]++].out = idx;
    //        m_incidentHalfEdges[currentPosition[m_indices[next]]++].in = idx;
    //    }
    //}
}

Mesh::Mesh(Mg::Trade::MeshData const& meshData) { setFromData(meshData); }

void Mesh::requireGaussianCurvature() {
    if(gaussianCurvature.size() != m_vertexCount) {
        requireAngles();
        arrayReserve(gaussianCurvature, m_vertexCount);
        for(Vertex vertex : vertices()) {
            Rad angleSum{0};
            for(Corner corner : vertex.corners())
                angleSum += corner.angle();
            gaussianCurvature[vertex] = 2.f*Mg::Math::Constants<float>::pi() - static_cast<float>(angleSum);
        }
    }
}

void Mesh::requireAngles() {
    if(angle.size() != m_cornerCount) {
        arrayResize(angle, m_cornerCount);
        for(Corner corner : corners()) {
            HalfEdge side1 = corner.side1();
            HalfEdge side2 = corner.side2();
            printf("half edge (%zu, %zu)\n", side1.tail().idx, side1.tip().idx);
            printf("half edge (%zu, %zu)\n", side2.tail().idx, side2.tip().idx);

            angle[corner] = Rad{Mg::Math::angle(side1.direction().normalized(), side2.direction().normalized())};
            Debug{} << Deg{angle[corner]};

            printf("\n");
        }
        Vertex v{4, this};
        for(Corner c : v.corners()) {
            printf("corner in face %zu has %f degree\n", c.face().idx, static_cast<float>(Deg{c.angle()}));
        }
    }
}

void Mesh::requireFaceAreas() {
    if(faceArea.size() != m_faceCount) {
        arrayResize(faceArea, m_faceCount);
        for(Face face : faces()) {
            Corner corner = face.halfEdge().corner();
            Vector3 side1 = corner.side1().direction();
            Vector3 side2 = corner.side2().direction();
            faceArea[face] = Mg::Math::cross(side1, side2).length()*0.5f;
        }
    }
}

void Mesh::requireEdgeLengths() {
    if(edgeLength.size() != m_edgeCount) {
        arrayResize(edgeLength, m_faceCount);
        for(Edge edge : edges()) {
            edgeLength[edge] = (edge.vertex1().position() - edge.vertex2().position()).length();
        }
    }
}

StridedArrayView1D<Vector3> Mesh::positions() {
    return stridedArrayView(m_attributes).slice(&Implementation::Attributes::position);
}

StridedArrayView1D<Vector3> Mesh::normals() {
    return stridedArrayView(m_attributes).slice(&Implementation::Attributes::normal);
}

StridedArrayView1D<Float> Mesh::scalars() {
    auto uvs = stridedArrayView(m_attributes).slice(&Implementation::Attributes::uv);
    return arrayCast<2, float>(uvs).slice<1>();
}

StridedArrayView1D<Color4> Mesh::colors() {
    return stridedArrayView(m_attributes).slice(&Implementation::Attributes::color);
}

//bool Mesh::isManifold() {
//    for (Edge e : edges()) {
//        if (!e.isManifold()) return false;
//    }
//    for (Vertex v : vertices()) {
//        if (!v.isManifold()) return false;
//    }
//    return true;
//}


VertexSet Mesh::vertices() { return {{0, this}, {m_vertexCount, this}}; }

FaceSet Mesh::faces() { return {{0, this}, {m_faceCount, this}}; }

EdgeSet Mesh::edges() { return {{0, this}, {m_edgeCount, this}}; }

HalfEdgeSet Mesh::halfEdges() { return {{0, this}, {m_halfEdgeCount, this}}; }

CornerSet Mesh::corners() {
    CornerIterator b{0, this};
    if(!b.isCorner()) ++b;
    return {b, {m_halfEdgeCount, this}};
}

void Mesh::uploadVertexBuffer(Mg::GL::Buffer& buffer) const { buffer.setData(m_attributes); }

void Mesh::uploadIndexBuffer(Mg::GL::Buffer& buffer) const { buffer.setData(m_indices); }

size_t Mesh::faceCount() const { return m_faceCount; }

size_t Mesh::vertexCount() const { return m_vertexCount; }

size_t Mesh::edgeCount() const { return m_edgeCount; }

size_t Mesh::indexCount() const { return m_indices.size(); }

size_t Mesh::halfEdgeCount() const { return m_halfEdgeCount; }

}
