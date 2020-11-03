//
// Created by janos on 9/5/20.
//

#include "MeshElements.h"
#include "Mesh.h"

#include <Corrade/Containers/StridedArrayView.h>
#include <Corrade/Containers/StaticArray.h>

#include <Magnum/Math/Vector3.h>
#include <Magnum/Math/Range.h>
#include <Magnum/Math/FunctionsBatch.h>

#include <cstdio>

namespace Phasefield {

/* ===== Half-Edge Implementation ===== */

HalfEdge HalfEdge::next() const { return HalfEdge{mesh->m_halfEdges[idx].next, mesh}; }

HalfEdge HalfEdge::twin() const { return HalfEdge{idx ^ 1, mesh}; }

Vertex HalfEdge::tip() const { return next().tail(); }

Vertex HalfEdge::tail() const { return Vertex{mesh->m_halfEdges[idx].tail, const_cast<Mesh*>(mesh)}; }

Face HalfEdge::face() const { return Face{mesh->m_halfEdges[idx].face, mesh}; }

Edge HalfEdge::edge() const { return Edge{idx/2, mesh}; }

Vector3 HalfEdge::direction() const { return tip().position() - tail().position(); }

HalfEdge HalfEdge::nextOutgoingHalfEdge() const { return twin().next(); }

HalfEdge HalfEdge::nextIncomingHalfEdge() const { return next().twin(); }

Corner HalfEdge::corner() const { return Corner{idx, mesh}; }

bool HalfEdge::isInterior() const { return !!face(); }

bool HalfEdge::onBoundaryLoop() const { return !face(); }

Vector3d HalfEdge::asVector() const { return Vector3d(tip().position() - tail().position()); }

Vector3d const& HalfEdge::gradient() const { return mesh->gradient[*this]; }

Debug& operator<<(Debug& debug, HalfEdge const& he) {
    char buffer[100];
    sprintf(buffer, "Half-Edge (%zu, %zu)", he.tail().idx, he.tip().idx);
    debug << buffer;
    return debug;
}

/* ===== Vertex ===== */

double Vertex::gaussianCurvature() const { return mesh->gaussianCurvature[*this]; }

Vector3 const& Vertex::normal() const { return mesh->normals()[idx]; }

Vector3 const& Vertex::position() const { return mesh->positions()[idx]; }

Float Vertex::scalar() const { return mesh->scalars()[idx]; }

HalfEdge Vertex::halfEdge() const { return HalfEdge{mesh->m_vertexHalfEdge[idx], mesh}; }

VertexCornerRange Vertex::corners() const {
    HalfEdge he = halfEdge().twin();
    if(he.onBoundaryLoop()) he = he.nextIncomingHalfEdge();
    return {{.he = he}, {.he = he}};
}

VertexVertexRange Vertex::adjacentVertices() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

VertexIncomingHalfEdgeRange Vertex::incomingHalfEdges() const { return {{.he = halfEdge().twin()}, {.he = halfEdge().twin()}}; }

VertexOutgoingHalfEdgeRange Vertex::outgoingHalfEdges() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

size_t Vertex::degree() const { return mesh->degree[*this]; }

VertexFaceRange Vertex::faces() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

VertexEdgeRange Vertex::edges() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

bool Vertex::onBoundary() const { return mesh->isOnBoundary[*this]; }

size_t Vertex::computeDegree() const {
    size_t count = 0;
    for(Vertex v : adjacentVertices())
        count++;
    return count;
}

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

double Face::area() const { return mesh->faceArea[*this]; }

Vector3d Face::normal() const { return mesh->faceNormal[*this]; }

FaceVertexRange Face::vertices() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

FaceHalfEdgeRange Face::halfEdges() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

FaceDualEdgeRange Face::dualEdges() const {
    HalfEdge he = halfEdge();
    do { he = he.next(); } while(!he.edge().hasDualEdge());
    return {{.he = he}, {.he = he}};
}

double Face::diameter() const { return mesh->faceDiameter[*this]; }

Range3D Face::bb() const {
    auto pos = positions();
    return Math::minmax({pos[0], pos[1], pos[2]});
}

StaticArray<3, Vector3> Face::positions() const {
    StaticArray<3, Vector3> pos;
    size_t i = 0;
    for(Vertex v : vertices())
        pos[i++] = v.position();
    return pos;
}

Vector3 Face::computeNormal() const {
    size_t i = 0;
    Vector3 positions[3];
    for(Vertex v : vertices()) positions[i++] = v.position();

    Vector3 side1 = positions[2] - positions[0];
    Vector3 side2 = positions[1] - positions[0];
    Vector3 normal = Math::cross(side1, side2);

    return normal;
}

bool Face::isTriangle() const {
    HalfEdge he = halfEdge();
    HalfEdge it = he;
    for(size_t i = 0; i < 3; ++i)
        it = it.next();
    return it == he;
}


Debug& operator<<(Debug& debug, Face const& f) {
    char buffer[100];
    Vertex vs[3];
    HalfEdge he = f.halfEdge();
    for(Vertex& v : vs) {
        v = he.tail();
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

bool Edge::hasDualEdge() const { return !!(halfEdge().face()) && !!(halfEdge().twin().face()); }

DualEdge Edge::dualEdge() const {
    CORRADE_ASSERT(hasDualEdge(), "edge does not have a dual edge", {});
    return {idx, mesh};
}

Vertex Edge::otherVertex(Vertex v) const {
    if(v == vertex1()) return vertex2();
    return vertex1();
}

double Edge::length() const { return halfEdge().asVector().length(); }

bool Edge::onBoundaryLoop() const { return halfEdge().onBoundaryLoop(); }

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

Radd Corner::angle() const { return mesh->angle[*this]; }

HalfEdge Corner::halfEdge() const { return HalfEdge{idx, mesh}; }

Rad Corner::computeAngle() const { return Mg::Math::angle(side1().direction().normalized(), side2().direction().normalized()); }

Debug& operator<<(Debug& debug, Corner const& c) {
    char buffer[100];
    HalfEdge s1 = c.side1();
    HalfEdge s2 = c.side2();
    sprintf(buffer, "Corner {{%zu, %zu}, {%zu, %zu}}", s1.tail().idx, s1.tip().idx, s2.tail().idx, s2.tip().idx);
    debug << buffer;
    return debug;
}

/* ===== Dual Edge ===== */

Face DualEdge::face1() const { return edge().halfEdge().face(); }

Face DualEdge::face2() const { return edge().halfEdge().twin().face(); }

Edge DualEdge::edge() const { return {idx, mesh}; }

Face DualEdge::otherFace(Face face) const {
    if(face == face1()) return face2();
    return face1();
}

Debug& operator<<(Debug& debug, DualEdge const& de) {
    char buffer[100];
    Face f1 = de.face1();
    Face f2 = de.face2();
    sprintf(buffer, "DualEdge {%zu, %zu}", f1.idx, f2.idx);
    debug << buffer;
    return debug;
}

}
