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
};

struct PhasefieldTree {
    Cr::Containers::Array<double> phasefieldData;
    Cr::Containers::Array<PhasefieldNode> nodes;
    uint32_t numLeafs = 0;
    uint32_t phasefieldSize = 0;
    uint32_t depth = 0;
    uint32_t firstLeafIdx = 0;

    PhasefieldNode& root() { return nodes.front(); }
    View2D<Mg::Double> phasefields() { return {phasefieldData, {nodes.size(), phasefieldSize}};}
    Range<PhasefieldNode*> leafs() { return {nodes.begin() + firstLeafIdx, nodes.end()}; }
    bool isLeftChild(int32_t idx) {
        auto parent = nodes[idx].parent;
        auto leftChildOfParent = nodes[parent].leftChild;
        return idx == leftChildOfParent;
    }
};
