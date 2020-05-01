//
// Created by janos on 2/23/20.
//

#pragma once

#include <Corrade/Containers/Array.h>

struct Node{
    int rank;
    int parent;
};

struct UnionFind : public Corrade::Containers::Array<Node> {

    explicit UnionFind(int n);

    [[nodiscard]] int& parent(int x);
    [[nodiscard]] int& rank(int x);

    int find(int x);

    void unite(int x, int y);

};

