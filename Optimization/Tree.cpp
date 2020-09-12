//
// Created by janos on 7/27/20.
//

#include "Tree.h"
#include "../Mesh/Mesh.h"

#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Utility/MurmurHash2.h>

#include <Magnum/Math/Vector4.h>
#include <Magnum/Magnum.h>

#include <cstdio>
#include <unordered_map>

namespace Phasefield {

namespace {

struct Edge {
    Edge(UnsignedInt v1_, UnsignedInt v2_) : v1(std::min(v1_, v2_)), v2(std::max(v1_, v2_)) {}

    UnsignedInt v1, v2;

    auto operator<=>(const Edge&) const = default;
};

struct EdgeHash {
    std::size_t operator()(Edge const& e) const {
        return *((std::size_t*) (Cr::Utility::MurmurHash2{}((char*) &e, sizeof(size_t)).byteArray()));
    }
};

}

Tree::Tree(Mesh& m) : mesh(&m), nodes(1), numLeafs(1) {}

//void Tree::subdivide(Containers::Array<UnsignedInt>& indices, Containers::Array<Vector3d>& vertices) {
//
//    std::unordered_map<Edge, int, EdgeHash> edgeMap(indices.size());
//
//    const auto pfs = phasefields();
//    Containers::Array<Containers::Array<double>> interpolated(phasefieldCount());
//
//    const auto vertexOffset = vertices.size();
//    int edgeCount = 0;
//
//    const std::size_t indexCount = indices.size();
//
//    for(std::size_t i = 0; i < indexCount; i += 3) {
//        UnsignedInt newVertices[3];
//        for(int j = 0; j < 3; ++j) {
//            const auto v1 = indices[i + j], v2 = indices[i + (j + 1)%3];
//            const auto[it, inserted] = edgeMap.try_emplace(Edge{v1, v2}, vertexOffset + edgeCount);
//            newVertices[j] = it->second;
//            if(inserted) {
//                Containers::arrayAppend(vertices, 0.5*(vertices[v1] + vertices[v2]));
//                for(std::size_t k = 0; k < phasefieldCount(); ++k) {
//                    Containers::arrayAppend(interpolated[k], 0.5*(pfs[j][v1] + pfs[j][v2]));
//                }
//                ++edgeCount;
//            }
//        }
//
//        /*
//            Add three new faces (0, 1, 3) and update original (2)
//
//                          orig 0
//                          /   \
//                         /  0  \
//                        /       \
//                    new 0 ----- new 2
//                    /   \       /  \
//                   /  1  \  2  / 3  \
//                  /       \   /      \
//             orig 1 ----- new 1 ---- orig 2
//        */
//        Containers::arrayAppend(indices, {
//                indices[i], newVertices[0], newVertices[2],
//                newVertices[0], indices[i + 1], newVertices[1],
//                newVertices[2], newVertices[1], indices[i + 2]});
//
//        for(std::size_t j = 0; j != 3; ++j)
//            indices[i + j] = newVertices[j];
//    }
//
//    /* flatten everything */
//    Containers::Array<double> data(phasefieldCount()*vertices.size());
//    Containers::StridedArrayView2D<double> phasefieldsNew(phasefieldData, {phasefieldCount(), vertices.size()});
//    for(std::size_t i = 0; i < phasefieldCount(); ++i) {
//        Utility::copy(pfs[i], phasefieldsNew[i].prefix(phasefieldSize));
//        Utility::copy(interpolated[i], phasefieldsNew[i].slice(phasefieldSize, vertices.size()));
//    }
//
//    /* mutate tree structure */
//    phasefieldData = std::move(data);
//    Containers::arrayResize(tempsData, phasefieldData.size());
//    phasefieldSize = vertices.size();
//}


void Tree::update() {
    size_t n = mesh->vertexCount();
    Array<double> data(DirectInit, n*nodeCount(), 0.5);
    StridedArrayView2D<double> view{data, {nodes.size(), n}};
    auto oldView = phasefields();
    for(std::size_t i = 0; i < nodes.size(); ++i) {
        size_t min = Math::min(n, m_vertexCount);
        Cr::Utility::copy(oldView[i].prefix(min), view[i].prefix(min));
    }

    m_vertexCount = n;
    arrayResize(tempsData, data.size());
    phasefieldData = std::move(data);
}

bool Tree::isLeftChild(Node const& node) const {
    auto leftChildOfParent = nodes[node.parent].leftChild;
    return node.idx == nodes[leftChildOfParent].idx;
}

Containers::StridedArrayView2D<Double> Tree::phasefields() {
    return {phasefieldData, {nodes.size(), m_vertexCount}};
}

Containers::StridedArrayView2D<Double> Tree::temporaryData() {
    return {tempsData, {nodes.size(), m_vertexCount}};
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "misc-no-recursion"

void Tree::remove(Node& nodeToDelete) {

    CORRADE_ASSERT(&nodeToDelete != &root(), "Cannot remove root node",);

    Containers::Array<double> newData;
    Containers::Array<Node> newNodes;

    size_t nodeIdx = 0;
    numLeafs = 0;
    depth = 0;

    auto visitor = YCombinator{
            [&](auto&& visitor, Node& node, size_t d = 0) -> size_t {
                if(&node == &nodeToDelete)
                    return Invalid;

                Int idx = nodeIdx++;

                Containers::ArrayView<const double> view{phasefields()[node.idx].asContiguous()};
                Containers::arrayAppend(newData, view);

                Node newNode = node;
                static_assert(size_t(-1) == Invalid);
                newNode.parent = idx - 1; /* For the root node this is -1 which is the same as Invalid */
                newNode.idx = idx;

                depth = Math::max(d, depth);

                if(!node.isLeaf()) {
                    if(node.leftChild != Invalid)
                        newNode.leftChild = visitor(nodes[node.leftChild], d + 1);
                    if(node.rightChild != Invalid)
                        newNode.rightChild = visitor(nodes[node.rightChild], d + 1);
                } else {
                    ++numLeafs;
                }

                Containers::arrayAppend(newNodes, newNode);

                return idx;
            }
    };

    visitor(root());

    Containers::arrayResize(tempsData, newData.size());
    phasefieldData = std::move(newData);

}

#pragma clang diagnostic pop

void Tree::addChild(Node& node, Child childType) {
    bool wasLeaf = node.isLeaf();

    arrayResize(phasefieldData, NoInit, (nodeCount() + 1)*vertexCount());
    arrayResize(tempsData, NoInit, (nodeCount() + 1)*vertexCount());

    size_t idx = levelStartIndex(node.depth + 1);
    double* data = phasefieldData.data();
    size_t sizeToMove = (nodes.size() - idx)*m_vertexCount*sizeof(double);
    if(sizeToMove) {
        memmove(data + (idx + 1)*m_vertexCount, data + idx*m_vertexCount, sizeToMove);
    }

    for(double& x : phasefields()[idx]) x = 0.;

    if(childType == Child::Right)
        node.rightChild = idx;
    else if(childType == Child::Left)
        node.leftChild = idx;

    Node child{.parent = node.idx, .idx = idx, .depth = node.depth + 1};
    arrayAppend(nodes, child);

    if(!wasLeaf)
        ++numLeafs;
    if(node.depth + 1 > depth)
        depth = node.depth + 1;

}

ArrayView<double> Tree::level(size_t d) {
    size_t begin = levelStartIndex(d);
    size_t end = levelStartIndex(d + 1);
    return phasefieldData.slice(begin*m_vertexCount, end*m_vertexCount);
}

void Tree::serialize(Array<char>& data) const {
    size_t dataSize = phasefieldData.size();
    arrayAppend(data, {(char*)&dataSize, sizeof(size_t)});
    arrayAppend(data, arrayCast<char>(phasefieldData));
    size_t nodesSize = nodes.size();
    arrayAppend(data, {(char*)&nodesSize, sizeof(size_t)});
    arrayAppend(data, arrayCast<char>(nodes));

    arrayAppend(data, {(char*)&numLeafs, sizeof(size_t)});
    arrayAppend(data, {(char*)&depth, sizeof(size_t)});
    arrayAppend(data, {(char*)&m_vertexCount, sizeof(size_t)});
}

Tree Tree::deserialize(Array<char> const& data, Mesh& m) {
    Tree t{m};
    char const* pc = data;
    size_t dataSize = deserializeTrivial<size_t>(pc);
    arrayResize(t.phasefieldData, NoInit, dataSize);
    arrayResize(t.tempsData, NoInit, dataSize);
    memcpy(t.phasefieldData.data(), pc, sizeof(double)*dataSize);
    pc += sizeof(double)*dataSize;

    size_t nodesSize = deserializeTrivial<size_t>(pc);
    arrayResize(t.nodes, NoInit, nodesSize);
    memcpy(t.nodes.data(), pc, sizeof(Node)*nodesSize);
    pc += sizeof(Node)*nodesSize;

    t.numLeafs = deserializeTrivial<size_t>(pc);
    t.depth = deserializeTrivial<size_t>(pc);
    t.m_vertexCount = deserializeTrivial<size_t>(pc);

    return t;
}

}