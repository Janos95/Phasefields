//
// Created by janos on 08.06.20.
//

#pragma once

#include "Enums.h"
#include "YCombinator.h"

#include <Corrade/Containers/Array.h>
#include <Magnum/Math/Functions.h>

namespace Phasefield {

namespace Cr = Corrade;

template<class Iter>
struct Range {
    Iter b, e;
    [[nodiscard]] auto begin() const { return b; }
    [[nodiscard]] auto end() const { return e; }
};

struct Node {
    static constexpr Mg::Int None = -1;

    Mg::Int leftChild = None;
    Mg::Int rightChild = None;
    Mg::Int parent = None;
    Mg::UnsignedInt idx;
    Mg::UnsignedInt depth;

    bool isLeaf() const {
        return leftChild == None && rightChild == None;
    }
};

struct Tree {

    Cr::Containers::Array<double> phasefieldData;
    Cr::Containers::Array<double> tempsData;

    Cr::Containers::Array<Node> nodes;
    uint32_t numLeafs = 0;
    uint32_t phasefieldSize = 0;
    uint32_t depth = 0;

    Node& root() { return nodes.front(); }

    Containers::StridedArrayView2D<Mg::Double> phasefields();

    Containers::StridedArrayView2D<Mg::Double> temps();

    Mg::UnsignedInt phasefieldCount() const {
        return nodes.size();
    }

    void resize(Mg::UnsignedInt n);

    void subdivide(Containers::Array<Mg::UnsignedInt>& indices, Containers::Array<Mg::Vector3d>& vertices);

    void remove(Node& node);

    void addLeftChild(Node& node);

    void addRightChild(Node& node);

    std::size_t phasefieldCount(Node& node){
        return nodes.size();
    }


    template<class F>
    void traverse(F&& f, Node* node = nullptr){
        auto visitor = YCombinator{
            [&](auto&& visitor, Node& node) -> void {
                f(node);
                if(!node.isLeaf()) { /* a leaf still has two children for computing derivatives */
                    if(node.leftChild != Node::None)
                        visitor(nodes[node.leftChild]);
                    if(node.rightChild != Node::None)
                        visitor(nodes[node.rightChild]);
                }
            }
        };
        visitor(node ? *node : root());
    }

    bool isLeftChild(Node const& node);
};

}