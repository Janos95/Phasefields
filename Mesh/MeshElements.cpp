//
// Created by janos on 9/5/20.
//

#include "MeshElements.h"
#include "Mesh.h"

#include <Corrade/Containers/StridedArrayView.h>
#include <Magnum/Math/Vector3.h>

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

bool HalfEdge::isInterior() const { return twin().face().isValid() && face().isValid(); }

bool HalfEdge::onBoundaryLoop() const { return !face().isValid(); };

/* ===== Vertex ===== */

Float Vertex::gaussianCurvature() const { return mesh->gaussianCurvature[*this]; }

Vector3& Vertex::normal() const { return mesh->normals()[idx]; }

Vector3& Vertex::position() const { return mesh->positions()[idx]; }

Float& Vertex::scalar() const { return mesh->scalars()[idx]; }

HalfEdge Vertex::halfEdge() const { return HalfEdge{mesh->m_vertexHalfEdge[idx], mesh}; }

VertexCornerRange Vertex::corners() const { return {{.he = halfEdge().twin()}, {.he = halfEdge().twin()}}; }

VertexAdjacentVertexRange Vertex::adjacentVertices() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

IncomingHalfEdgeRange Vertex::incomingHalfedges() const { return {{.he = halfEdge().twin()}, {.he = halfEdge().twin()}}; }

OutgoingHalfEdgeRange Vertex::outgoingHalfedges() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }


/* ===== Face ===== */

HalfEdge Face::halfEdge() const { return HalfEdge{mesh->m_faceHalfEdge[idx], mesh}; }

CornersOfFaceRange Face::corners() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

FaceEdgeRange Face::edges() const { return {{.he = halfEdge()}, {.he = halfEdge()}}; }

bool Face::isValid() const { return idx != Invalid; }

/* ===== Edge ===== */

HalfEdge Edge::halfEdge() { return HalfEdge{2*idx, mesh}; }

Vertex Edge::vertex1() { return halfEdge().tail(); }

Vertex Edge::vertex2() { return halfEdge().tip(); }

/* ===== Corner ===== */

HalfEdge Corner::side1() const { return halfEdge().twin(); }

HalfEdge Corner::side2() const { return halfEdge().next(); }

Vertex Corner::vertex() const { return halfEdge().tip(); }

Face Corner::face() const { return halfEdge().face(); }

Rad Corner::angle() const { return mesh->angle[*this]; }

HalfEdge Corner::halfEdge() const { return HalfEdge{idx, mesh}; }

}
