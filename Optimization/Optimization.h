//
// Created by janos on 7/1/20.
//

#pragma once

//forward declarations

namespace Phasefield {

struct Functional;
struct LossFunction;
struct SparseMatrix;

struct Tree;
struct Node;

template<int IteratorType> struct NodeIterator;
using LeafIterator = NodeIterator<0>;
using InternalNodeIterator = NodeIterator<1>;
using HorizontalNodeIterator = NodeIterator<2>;

namespace Solver {

struct Problem;
struct RecursiveProblem;

struct IterationSummary;

struct Options;
struct Summary;

}

}