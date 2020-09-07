//
// Created by janos on 9/7/20.
//

#pragma once

namespace Phasefield {

class Mesh;

struct Vertex;
struct Edge;
struct HalfEdge;
struct Face;
struct Corner;

template<class>
struct ElementIterator;

struct VertexCornerIterator;
struct VertexAdjacentVertexIterator;
struct IncomingHalfEdgeIterator;
struct OutgoingHalfEdgeIterator;
struct FaceEdgeIterator;
struct CornersOfFaceIterator;

using VertexIterator = ElementIterator<Vertex>;
using FaceIterator = ElementIterator<Face>;
using EdgeIterator = ElementIterator<Edge>;
using HalfEdgeIterator = ElementIterator<HalfEdge>;
using CornerIterator = ElementIterator<Corner>;

template<class> struct Range;

using VertexCornerRange = Range<VertexCornerIterator>;
using VertexAdjacentVertexRange = Range<VertexAdjacentVertexIterator>;
using IncomingHalfEdgeRange = Range<IncomingHalfEdgeIterator>;
using OutgoingHalfEdgeRange = Range<OutgoingHalfEdgeIterator>;
using FaceEdgeRange = Range<FaceEdgeIterator>;
using CornersOfFaceRange = Range<CornersOfFaceIterator>;

using VertexSet = Range<VertexIterator>;
using FaceSet = Range<FaceIterator>;
using EdgeSet = Range<EdgeIterator>;
using HalfEdgeSet = Range<HalfEdgeIterator>;
using CornerSet = Range<CornerIterator>;

}