//
// Created by janos on 2/23/20.
//

#include "UnionFind.h"

namespace Phasefield {

using Node = Implementation::Node;

UnionFind::UnionFind(size_t n) : Array<Implementation::Node>(NoInit, n) {
    for(size_t i = 0; i < n; ++i) {
        (*this)[i] = Node{1, i};
    }
}

size_t& UnionFind::parent(size_t x) {
    return (*this)[x].parent;
}

size_t& UnionFind::rank(size_t x) {
    return (*this)[x].rank;
}


size_t UnionFind::find(size_t x) {
    size_t root = x;
    while(parent(root) != root) {
        root = parent(root);
    }

    while(parent(x) != root) {
        auto p = parent(x);
        (*this)[x].parent = root;
        x = p;
    }

    return root;
}

void UnionFind::unite(size_t x, size_t y) {
    auto xRoot = find(x);
    auto yRoot = find(y);

    if(xRoot == yRoot)
        return;

    if(rank(xRoot) < rank(yRoot))
        std::swap(xRoot, yRoot);

    parent(yRoot) = xRoot; //
    if(rank(xRoot) == rank(yRoot))
        ++rank(xRoot);
}

}