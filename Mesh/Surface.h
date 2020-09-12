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

template<class> struct ElementIterator;
template<class> struct FaceCirculationIterator;

struct VertexCornerIterator;
struct VertexAdjacentVertexIterator;
struct IncomingHalfEdgeIterator;
struct OutgoingHalfEdgeIterator;

using FaceCornerIterator = FaceCirculationIterator<Corner>;
using FaceAdjacentFaceIterator = FaceCirculationIterator<Face>;
using FaceVertexIterator = FaceCirculationIterator<Vertex>;
using FaceEdgeIterator = FaceCirculationIterator<Edge>;
using FaceHalfEdgeIterator = FaceCirculationIterator<HalfEdge>;

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
using FaceCornerRange = Range<FaceCornerIterator>;
using FaceAdjacentFaceRange = Range<FaceAdjacentFaceIterator>;
using FaceVertexRange = Range<FaceVertexIterator>;
using FaceHalfEdgeRange = Range<FaceHalfEdgeIterator>;

using VertexSet = Range<VertexIterator>;
using FaceSet = Range<FaceIterator>;
using EdgeSet = Range<EdgeIterator>;
using HalfEdgeSet = Range<HalfEdgeIterator>;
using CornerSet = Range<CornerIterator>;

}