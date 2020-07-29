//
// Created by janos on 08.06.20.
//

#pragma once

#include "types.hpp"
#include "y_combinator.hpp"

#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/StridedArrayView.h>
#include <Corrade/Utility/Algorithms.h>
#include <Magnum/Math/Functions.h>

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
        return {tempsData, {nodes.size(), phasefieldSize}};
    }

    Mg::UnsignedInt phasefieldCount() const{
        return nodes.size();
    }

    void resize(Mg::UnsignedInt n){
        Containers::Array<double> data(n * phasefieldSize);
        Containers::StridedArrayView2D<double> view{data, {nodes.size(), n}};
        auto oldView = phasefields();
        for(std::size_t i = 0; i < nodes.size(); ++i)
            Cr::Utility::copy(oldView[i].slice(0, Mg::Math::min(n, phasefieldSize)), view[i]);
        phasefieldSize = n;

        Containers::arrayResize(tempsData, data.size());
        phasefieldData = std::move(data);
    }

    void subdivide(Containers::Array<Mg::UnsignedInt>& indices, Containers::Array<Mg::Vector3d>& vertices);

    template<class F>
    void traverse(F&& f){
        auto visitor = YCombinator{
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
