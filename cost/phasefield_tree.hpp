//
// Created by janos on 08.06.20.
//

#pragma once

#include "types.hpp"

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/StridedArrayView.h>

namespace Cr = Corrade;

template<class Iter>
struct Range {
    Iter b, e;
    [[nodiscard]] auto begin() const { return b; }
    [[nodiscard]] auto end() const { return e; }
};

struct PhasefieldNode {
    static constexpr int32_t None = -1;

    int32_t leftChild = None;
    int32_t rightChild = None;
    int32_t parent = None;
    uint32_t idx;
    uint32_t depth;

    bool isLeaf() const {
        return leftChild == None && rightChild == None;
    }
};

struct PhasefieldTree {
    Cr::Containers::Array<double> phasefieldData;
    Cr::Containers::Array<double> tempsData;

    Cr::Containers::Array<PhasefieldNode> nodes;
    uint32_t numLeafs = 0;
    uint32_t phasefieldSize = 0;
    uint32_t depth = 0;

    PhasefieldNode& root() { return nodes.front(); }

    Containers::StridedArrayView2D<Mg::Double> phasefields() {
        return {phasefieldData, {nodes.size(), phasefieldSize}};
    }

    Containers::StridedArrayView2D<Mg::Double> temps() {
        return {tempsData, {tempsData.size(), phasefieldSize}};
    }

    template<class F>
    void traverse(F&& f){
        auto visitor = YCombinator {
            [&](auto&& visitor, PhasefieldNode& node) -> void {
                f(node);
                if(!node.isLeaf()) { /* a leaf still has two children for computing derivatives */
                    if(node.leftChild != PhasefieldNode::None)
                        visitor(nodes[node.leftChild]);
                    if(node.rightChild != PhasefieldNode::None)
                        visitor(nodes[node.rightChild]);
                }
            }
        };
        visitor(root());
    }

    bool isLeftChild(PhasefieldNode const& node) {
        auto leftChildOfParent = nodes[node.parent].leftChild;
        return node.idx == nodes[leftChildOfParent].idx;
    }
};
