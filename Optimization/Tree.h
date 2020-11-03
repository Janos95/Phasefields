//
// Created by janos on 08.06.20.
//

#pragma once

#include "Enums.h"
#include "Range.h"
#include "Types.h"
#include "Mesh.h"
#include "Serialize.h"
#include "Optimization.h"

#include <Corrade/Containers/Array.h>
#include <Magnum/Math/Functions.h>

namespace Phasefield {

namespace Cr = Corrade;

struct NodeData {
    size_t leftChild = Invalid;
    size_t rightChild = Invalid;
    size_t parent = Invalid;
    size_t depth = 0;
};

struct Node {
    size_t idx = Invalid;
    Tree* tree = nullptr;

    auto operator<=>(Node const& other) const { return idx <=> other.idx; }
    bool operator==(Node const& other) const = default;

    [[nodiscard]] bool isLeaf() const { return !hasLeftChild() && !hasRightChild(); }

    [[nodiscard]] bool isValid() const { return idx != Invalid; }

    [[nodiscard]] bool hasLeftChild() const { return leftChild().isValid(); }

    [[nodiscard]] bool hasRightChild() const { return rightChild().isValid(); }

    [[nodiscard]] bool isLeftChild() const { return parent().leftChild() == *this; }

    [[nodiscard]] bool isRightChild() const { return parent().rightChild() == *this; }

     void splitAndInitialize(Node* n);

    [[nodiscard]] double integrateWeight(Mesh& mesh) const;

    [[nodiscard]] VertexDataView<Double> phasefield() const;

    [[nodiscard]] VertexDataView<Double> temporary() const;

    [[nodiscard]] Node parent() const;

    Node addLeftChild(Node* n = nullptr);

    Node addRightChild(Node* n = nullptr);

    Node addChild(bool left, Node* n);

    void initializePhasefieldFromParent();

    [[nodiscard]] size_t depth() const;

    [[nodiscard]] Node leftChild() const;

    [[nodiscard]] Node rightChild() const;

    explicit operator bool() const { return idx != Invalid; }
};

Debug& operator<<(Debug& debug, Node const& n);

struct Tree {

    friend Node;

    explicit Tree(Mesh&);

    Mesh* mesh = nullptr;

    Array<double> phasefieldData;
    Array<double> tempsData;

    Array<NodeData> nodeData;

    size_t numLeafs = 0;
    size_t depth = 0;

    Node root() { return {0, this}; }

    ArrayView<double> level(size_t d);

    void update();

    /**
     * If the level does not exist, the number of nodes is returned.
     * @param level depth in the tree
     * @return idx of the first node on the requested level
     */
    size_t levelStartIndex(size_t level);

    size_t nodeCountOnLevel(size_t level) { return levelStartIndex(level + 1) - levelStartIndex(level); }

    //void subdivide(Containers::Array<Mg::UnsignedInt>& indices, Containers::Array<Vector3d>& vertices);

    //void remove(Node node);
    StridedArrayView2D<double> phasefields();

    Node insertNodeAtIndex(size_t idx);

    void reset();

    Range<HorizontalNodeIterator> nodesOnLevel(size_t l);

    Range<HorizontalNodeIterator> nodesBelowLevel(size_t l);

    Range<HorizontalNodeIterator> nodes();

    Range<InternalNodeIterator> internalNodes();

    Array<Node> ancestorsOfLevel(size_t l);

    Range<LeafIterator> leafs();

    [[nodiscard]] size_t nodeCount() const { return nodeData.size(); }

    [[nodiscard]] size_t vertexCount() const { return phasefieldData.size()/nodeCount(); }

    void computeWeightsOfAncestorsOfLevel(size_t l);

    void computeLeafWeights();

    void serialize(Array<char>& data) const;

    static Tree deserialize(ArrayView<const char> const& data, Mesh& m);
};

template<int IteratorType>
struct NodeIterator {
    Node node;

    Node operator *() const { return node; }
    NodeIterator& operator++() requires (IteratorType == 2){ node.idx++; return *this; }

    NodeIterator& operator++() requires (IteratorType == 0) {
        do { node.idx++; } while(node.idx < node.tree->nodeCount() && !node.isLeaf());
        return *this;
    }

    NodeIterator& operator++() requires (IteratorType == 1) {
        do { node.idx++; } while(node.idx < node.tree->nodeCount() && node.isLeaf());
        return *this;
    }

    bool operator !=(NodeIterator const& other) const { return node != other.node; }
};


}