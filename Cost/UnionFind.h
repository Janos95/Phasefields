//
// Created by janos on 2/23/20.
//

#pragma once

#include "Types.h"
#include <Corrade/Containers/Array.h>

namespace Phasefield {

namespace Implementation {
struct Node {
    size_t rank;
    size_t parent;
};

}

struct UnionFind : Array<Implementation::Node> {

    explicit UnionFind(size_t n);

    [[nodiscard]] size_t& parent(size_t x);

    [[nodiscard]] size_t& rank(size_t x);

    size_t find(size_t x);

    void unite(size_t x, size_t y);

};

}

