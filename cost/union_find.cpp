//
// Created by janos on 2/23/20.
//

#include "union_find.hpp"

using namespace Corrade;

UnionFind::UnionFind(int n): Containers::Array<Node>(Containers::NoInit, n)
{
    for (int i = 0; i < n; ++i) {
        (*this)[i] = Node{1, i};
    }
}

int& UnionFind::parent(int x) {
    return (*this)[x].parent;
}

int& UnionFind::rank(int x) {
    return (*this)[x].rank;
}


int UnionFind::find(int x) {
    auto root = x;
    while(parent(root) != root){
        root = parent(root);
    }

    while(parent(x) != root){
        auto p = parent(x);
        (*this)[x].parent = root;
        x = p;
    }

    return root;
}

void UnionFind::unite(int x, int y){
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

