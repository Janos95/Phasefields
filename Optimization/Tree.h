//
// Created by janos on 08.06.20.
//

#pragma once

#include "Enums.h"
#include "YCombinator.h"
#include "Range.h"
#include "Types.h"
#include "Mesh.h"

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

    bool isLeaf() const {
        return leftChild == Invalid && rightChild == Invalid;
    }
};

struct Tree {

    Mesh& mesh;

    Cr::Containers::Array<double> phasefieldData;
    Cr::Containers::Array<double> tempsData;

    Cr::Containers::Array<Node> nodes;
    size_t numLeafs = 0;
    size_t phasefieldSize = 0;
    size_t depth = 0;

    Node& root() { return nodes.front(); }

    Containers::StridedArrayView2D<Mg::Double> phasefields();

    Containers::StridedArrayView2D<Mg::Double> temps();

    Mg::UnsignedInt phasefieldCount() const {
        return nodes.size();
    }

    void update();

    //void subdivide(Containers::Array<Mg::UnsignedInt>& indices, Containers::Array<Vector3d>& vertices);

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
                    if(node.leftChild != Invalid)
                        visitor(nodes[node.leftChild]);
                    if(node.rightChild != Invalid)
                        visitor(nodes[node.rightChild]);
                }
            }
        };
        visitor(node ? *node : root());
    }

    bool isLeftChild(Node const& node);
};

}