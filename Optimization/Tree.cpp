//
// Created by janos on 7/27/20.
//

#include "Tree.h"

#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Utility/MurmurHash2.h>

#include <Magnum/Math/Vector4.h>
#include <Magnum/Magnum.h>

#include <cstdio>
#include <unordered_map>

namespace Phasefield {

using namespace Magnum;
using namespace Corrade;

namespace {

struct Edge {
    Edge(Mg::UnsignedInt v1_, Mg::UnsignedInt v2_) : v1(std::min(v1_, v2_)), v2(std::max(v1_, v2_)) {}

    Mg::UnsignedInt v1, v2;

    auto operator<=>(const Edge&) const = default;
};

struct EdgeHash {
    std::size_t operator()(Edge const& e) const {
        return *((std::size_t*) (Utility::MurmurHash2{}((char*) &e, sizeof(UnsignedLong)).byteArray()));
    }
};

}

void Tree::subdivide(Containers::Array<Mg::UnsignedInt>& indices, Containers::Array<Mg::Vector3d>& vertices) {

    std::unordered_map<Edge, int, EdgeHash> edgeMap(indices.size());

    const auto pfs = phasefields();
    Containers::Array<Containers::Array<double>> interpolated(phasefieldCount());

    const auto vertexOffset = vertices.size();
    int edgeCount = 0;

    const std::size_t indexCount = indices.size();

    for(std::size_t i = 0; i < indexCount; i += 3) {
        UnsignedInt newVertices[3];
        for(int j = 0; j < 3; ++j) {
            const auto v1 = indices[i + j], v2 = indices[i + (j + 1)%3];
            const auto[it, inserted] = edgeMap.try_emplace(Edge{v1, v2}, vertexOffset + edgeCount);
            newVertices[j] = it->second;
            if(inserted) {
                Containers::arrayAppend(vertices, 0.5*(vertices[v1] + vertices[v2]));
                for(std::size_t k = 0; k < phasefieldCount(); ++k) {
                    Containers::arrayAppend(interpolated[k], 0.5*(pfs[j][v1] + pfs[j][v2]));
                }
                ++edgeCount;
            }
        }

        /*
            Add three new faces (0, 1, 3) and update original (2)

                          orig 0
                          /   \
                         /  0  \
                        /       \
                    new 0 ----- new 2
                    /   \       /  \
                   /  1  \  2  / 3  \
                  /       \   /      \
             orig 1 ----- new 1 ---- orig 2
        */
        Containers::arrayAppend(indices, {
                indices[i], newVertices[0], newVertices[2],
                newVertices[0], indices[i + 1], newVertices[1],
                newVertices[2], newVertices[1], indices[i + 2]});

        for(std::size_t j = 0; j != 3; ++j)
            indices[i + j] = newVertices[j];
    }

    /* flatten everything */
    Containers::Array<double> data(phasefieldCount()*vertices.size());
    Containers::StridedArrayView2D<double> phasefieldsNew(phasefieldData, {phasefieldCount(), vertices.size()});
    for(std::size_t i = 0; i < phasefieldCount(); ++i) {
        Utility::copy(pfs[i], phasefieldsNew[i].prefix(phasefieldSize));
        Utility::copy(interpolated[i], phasefieldsNew[i].slice(phasefieldSize, vertices.size()));
    }

    /* mutate tree structure */
    phasefieldData = std::move(data);
    Containers::arrayResize(tempsData, phasefieldData.size());
    phasefieldSize = vertices.size();
}


void Tree::resize(Mg::UnsignedInt n) {
    Containers::Array<double> data(n*phasefieldSize);
    Containers::StridedArrayView2D<double> view{data, {nodes.size(), n}};
    auto oldView = phasefields();
    for(std::size_t i = 0; i < nodes.size(); ++i)
        Cr::Utility::copy(oldView[i].slice(0, Mg::Math::min(n, phasefieldSize)), view[i]);
    phasefieldSize = n;

    Containers::arrayResize(tempsData, data.size());
    phasefieldData = std::move(data);
}

bool Tree::isLeftChild(Node const& node) {
    auto leftChildOfParent = nodes[node.parent].leftChild;
    return node.idx == nodes[leftChildOfParent].idx;
}

Containers::StridedArrayView2D<Mg::Double> Tree::phasefields() {
    return {phasefieldData, {nodes.size(), phasefieldSize}};
}

Containers::StridedArrayView2D<Mg::Double> Tree::temps() {
    return {tempsData, {nodes.size(), phasefieldSize}};
}

#pragma clang diagnostic push
#pragma ide diagnostic ignored "misc-no-recursion"

void Tree::remove(Node& nodeToDelete) {

    CORRADE_ASSERT(&nodeToDelete != &root(), "Cannot remove root node",);


    Containers::Array<double> newData;
    Containers::Array<Node> newNodes;

    int nodeIdx = 0;
    numLeafs = 0;
    depth = 0;

    auto visitor = YCombinator{
            [&](auto&& visitor, Node& node, UnsignedInt d = 0) -> Int {
                if(&node == &nodeToDelete)
                    return Node::None;

                Int idx = nodeIdx++;

                Containers::ArrayView<const double> view{phasefields()[node.idx].asContiguous()};
                Containers::arrayAppend(newData, view);

                Node newNode = node;
                newNode.parent = idx - 1; /* For the root node this is -1 which is the same as PhasefieldNode::None */
                newNode.idx = idx;

                depth = Math::max(d, depth);

                if(!node.isLeaf()) {
                    if(node.leftChild != Node::None)
                        newNode.leftChild = visitor(nodes[node.leftChild], d + 1);
                    if(node.rightChild != Node::None)
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

}