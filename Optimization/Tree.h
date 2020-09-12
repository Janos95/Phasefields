//
// Created by janos on 08.06.20.
//

#pragma once

#include "Enums.h"
#include "YCombinator.h"
#include "Range.h"
#include "Types.h"
#include "Mesh.h"
#include "Serialize.h"

#include <Corrade/Containers/Array.h>
#include <Magnum/Math/Functions.h>

namespace Phasefield {

namespace Cr = Corrade;

struct Node {

    size_t leftChild = Invalid;
    size_t rightChild = Invalid;
    size_t parent = Invalid;
    size_t idx;
    size_t depth;

    [[nodiscard]] bool isLeaf() const {
        return leftChild == Invalid && rightChild == Invalid;
    }
};

struct Tree {

    enum class Child {
        Left,
        Right
    };

    explicit Tree(Mesh&);

    Mesh* mesh = nullptr;

    Array<double> phasefieldData;
    Array<double> tempsData;

    Array<Node> nodes;
    size_t numLeafs = 0;
    size_t depth = 0;
    size_t m_vertexCount = 0;

    Node& root() { return nodes.front(); }

    ArrayView<double> level(size_t d);

    StridedArrayView2D<Double> phasefields();

    StridedArrayView2D<Double> temporaryData();

    void update();

    /**
     * If the level does not exist, the number of nodes is returned.
     * @param level depth in the tree
     * @return idx of the first node on the requested level
     */
    size_t levelStartIndex(size_t level) {
        size_t nodeIdx = nodes.size();
        traverse([&](Node& node){
            if(node.depth == level) {
                if(nodeIdx == nodes.size())
                    nodeIdx = &node - nodes.begin();
                return false;
            } else return true;
        });
        return nodeIdx;
    }

    size_t nodeCountOnLevel(size_t level) {
        return levelStartIndex(level + 1) - levelStartIndex(level);
    }

    //void subdivide(Containers::Array<Mg::UnsignedInt>& indices, Containers::Array<Vector3d>& vertices);

    void remove(Node& node);

    void addChild(Node& node, Child child);

    [[nodiscard]] size_t nodeCount() const { return nodes.size(); }

    [[nodiscard]] size_t vertexCount() const { return m_vertexCount; }

    template<class F>
    void traverse(F&& f, Node* node = nullptr){
        auto visitor = YCombinator{
            [&](auto&& visitor, Node& node) -> void {
                if(!f(node)) return;
                if(!node.isLeaf()) {
                    if(node.leftChild != Invalid)
                        visitor(nodes[node.leftChild]);
                    if(node.rightChild != Invalid)
                        visitor(nodes[node.rightChild]);
                }
            }
        };
        visitor(node ? *node : root());
    }

    [[nodiscard]] bool isLeftChild(Node const& node) const;

    void serialize(Array<char>& data) const;

    static Tree deserialize(Array<char> const& data, Mesh& m);
};

}