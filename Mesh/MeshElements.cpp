//
// Created by janos on 9/5/20.
//

#include "MeshElements.h"
#include "Mesh.h"

#include <Corrade/Containers/StridedArrayView.h>
#include <Magnum/Math/Vector3.h>

#include <cstdio>

namespace Phasefield {

/* ===== Half-Edge Implementation ===== */

HalfEdge HalfEdge::next() const { return HalfEdge{mesh->m_halfEdges[idx].next, mesh}; }

HalfEdge HalfEdge::twin() const { return HalfEdge{idx ^ 1, mesh}; }

Vertex HalfEdge::tip() const { return next().tail(); }

Vertex HalfEdge::tail() const { return Vertex{mesh->m_halfEdges[idx].tail, mesh}; }

Face HalfEdge::face() const { return Face{mesh->m_halfEdges[idx].face, mesh}; }

Edge HalfEdge::edge() const { return Edge{idx/2, mesh}; }

Vector3 HalfEdge::direction() const { return tip().position() - tail().position(); }

HalfEdge HalfEdge::nextOutgoingHalfEdge() const { return twin().next(); }

HalfEdge HalfEdge::nextIncomingHalfEdge() const { return next().twin(); }

Corner HalfEdge::corner() const { return Corner{idx, mesh}; }

bool HalfEdge::isInterior() const { return face().isValid(); }

bool HalfEdge::onBoundaryLoop() const { return !face().isValid(); }

Vector3d HalfEdge::asVector() const { return Vector3d(tip().position() - tail().position()); }

Debug& operator<<(Debug& debug, HalfEdge const& he) {
    char buffer[100];
    sprintf(buffer, "Half-Edge (%zu, %zu)", he.tail().idx, he.tip().idx);
    debug << buffer;
    return debug;
}

/* ===== Vertex ===== */

Float Vertex::gaussianCurvature() const { return mesh->gaussianCurvature[*this]; }

Vector3& Vertex::normal() const { return mesh->normals()[idx]; }

Vector3& Vertex::position() const { return mesh->positions()[idx]; }

Float& Vertex::scalar() const { return mesh->scalars()[idx]; }

HalfEdge Vertex::halfEdge() const { return HalfEdge{mesh->m_vertexHalfEdge[idx], mesh}; }

VertexCornerRange Vertex::corners() const { return {{.he = halfEdge().twin()}, {.he = halfEdge().twin()}}; }

VertexAdjacentVertexRange Vertex::adjacentVertices() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

IncomingHalfEdgeRange Vertex::incomingHalfEdges() const { return {{.he = halfEdge().twin()}, {.he = halfEdge().twin()}}; }

OutgoingHalfEdgeRange Vertex::outgoingHalfEdges() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

Debug& operator<<(Debug& debug, Vertex const& v) {
    char buffer[100];
    sprintf(buffer, "Vertex %zu", v.idx);
    debug << buffer;
    return debug;
}

/* ===== Face ===== */

HalfEdge Face::halfEdge() const { return HalfEdge{mesh->m_faceHalfEdge[idx], mesh}; }

FaceCornerRange Face::corners() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

FaceEdgeRange Face::edges() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

bool Face::isValid() const { return idx != Invalid; }

double Face::area() const { return mesh->faceArea[*this]; }

Vector3d Face::normal() const { return mesh->faceNormal[*this]; }

FaceAdjacentFaceRange Face::adjacentFaces() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

FaceVertexRange Face::vertices() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

FaceHalfEdgeRange Face::halfEdges() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

Debug& operator<<(Debug& debug, Face const& f) {
    char buffer[100];
    Vertex vs[3];
    HalfEdge he = f.halfEdge();
    for(size_t i = 0; i < 3; ++i) {
        vs[i] = he.tail();
        he = he.next();
    }
    sprintf(buffer, "Face {%zu, %zu, %zu}", vs[0].idx, vs[1].idx, vs[2].idx);
    debug << buffer;
    return debug;
}

/* ===== Edge ===== */

HalfEdge Edge::halfEdge() const { return HalfEdge{2*idx, mesh}; }

Vertex Edge::vertex1() const { return halfEdge().tail(); }

Vertex Edge::vertex2() const { return halfEdge().tip(); }

Debug& operator<<(Debug& debug, Edge const& e) {
    char buffer[100];
    sprintf(buffer, "Edge {%zu, %zu}", e.vertex1().idx, e.vertex2().idx);
    debug << buffer;
    return debug;
}

/* ===== Corner ===== */

HalfEdge Corner::side1() const { return halfEdge().twin(); }

HalfEdge Corner::side2() const { return halfEdge().next(); }

Vertex Corner::vertex() const { return halfEdge().tip(); }

Face Corner::face() const { return halfEdge().face(); }

Rad Corner::angle() const { return mesh->angle[*this]; }

HalfEdge Corner::halfEdge() const { return HalfEdge{idx, mesh}; }

Debug& operator<<(Debug& debug, Corner const& c) {
    char buffer[100];
    HalfEdge s1 = c.side1();
    HalfEdge s2 = c.side2();
    sprintf(buffer, "Corner {{%zu, %zu}, {%zu, %zu}}", s1.tail().idx, s1.tip().idx, s2.tail().idx, s2.tip().idx);
    debug << buffer;
    return debug;
}

}
