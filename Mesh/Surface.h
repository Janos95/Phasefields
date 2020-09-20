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
struct DualEdge;

template<class> struct ElementIterator;
template<class> struct FaceCirculationIterator;

struct VertexCornerIterator;
struct VertexAdjacentVertexIterator;
struct IncomingHalfEdgeIterator;
struct OutgoingHalfEdgeIterator;

using FaceCornerIterator = FaceCirculationIterator<Corner>;
using FaceVertexIterator = FaceCirculationIterator<Vertex>;
using FaceEdgeIterator = FaceCirculationIterator<Edge>;
using FaceHalfEdgeIterator = FaceCirculationIterator<HalfEdge>;
using FaceDualEdgeIterator = FaceCirculationIterator<DualEdge>;

using VertexIterator = ElementIterator<Vertex>;
using FaceIterator = ElementIterator<Face>;
using EdgeIterator = ElementIterator<Edge>;
using HalfEdgeIterator = ElementIterator<HalfEdge>;
using CornerIterator = ElementIterator<Corner>;
using DualEdgeIterator = ElementIterator<DualEdge>;

template<class> struct Range;

using VertexCornerRange = Range<VertexCornerIterator>;
using VertexAdjacentVertexRange = Range<VertexAdjacentVertexIterator>;
using IncomingHalfEdgeRange = Range<IncomingHalfEdgeIterator>;
using OutgoingHalfEdgeRange = Range<OutgoingHalfEdgeIterator>;

using FaceEdgeRange = Range<FaceEdgeIterator>;
using FaceCornerRange = Range<FaceCornerIterator>;
using FaceVertexRange = Range<FaceVertexIterator>;
using FaceHalfEdgeRange = Range<FaceHalfEdgeIterator>;
using FaceDualEdgeRange = Range<FaceDualEdgeIterator>;

using VertexSet = Range<VertexIterator>;
using FaceSet = Range<FaceIterator>;
using EdgeSet = Range<EdgeIterator>;
using HalfEdgeSet = Range<HalfEdgeIterator>;
using CornerSet = Range<CornerIterator>;
using DualEdgeSet = Range<DualEdgeIterator>;


template<class E, class T> class MeshDataView;
template<class E, class T> class MeshData;

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

}