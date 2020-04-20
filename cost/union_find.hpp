//
// Created by janos on 2/23/20.
//

#pragma once

#include <Corrade/Containers/Array.h>

struct Node{
    int rank;
    int parent;
};

class UnionFind : public Corrade::Containers::Array<Node> {
public:

    explicit UnionFind(int n);

    int parent(int x) const;

    int find(int x);

    int unite(int x, int y);

    bool isValid(int k) const;

};

