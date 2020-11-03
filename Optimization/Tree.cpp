//
// Created by janos on 7/27/20.
//

#include "Tree.h"
#include "Mesh.h"
#include "C1Functions.h"
#include "Bfs.h"
#include "Bvh.h"

#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Utility/MurmurHash2.h>
#include <Corrade/Utility/FormatStl.h>

#include <Magnum/Math/Vector4.h>
#include <Magnum/Magnum.h>

#include <cstdio>
#include <unordered_map>
#include <queue>

namespace Phasefield {


Node Node::parent() const { return {tree->nodeData[idx].parent, tree}; }

Node Node::leftChild() const { return {tree->nodeData[idx].leftChild, tree}; }

Node Node::rightChild() const { return {tree->nodeData[idx].rightChild, tree}; }

size_t Node::depth() const { return tree->nodeData[idx].depth; }

VertexDataView<Double> Node::phasefield() const {
    size_t n = tree->vertexCount();
    return tree->phasefieldData.slice(idx*n, (idx + 1)*n);
}

VertexDataView<Double> Node::temporary() const {
    size_t n = tree->vertexCount();
    return tree->tempsData.slice(idx*n, (idx + 1)*n);
}

Node Node::addChild(bool left, Node* n) {

    bool wasLeaf = isLeaf();

    size_t i = tree->levelStartIndex(depth() + 1);
    for(Node node : tree->nodesOnLevel(depth() + 1)) {
        if(node.parent() > *this || (node.parent() == *this && left))
            break;
        i++;
    }

    Node child = tree->insertNodeAtIndex(i);

    /* hook up both pointers */
    tree->nodeData[child.idx].parent = idx;
    if(left) {
        tree->nodeData[idx].leftChild = child.idx;
    } else {
        tree->nodeData[idx].rightChild = child.idx;
    }
    tree->nodeData[child.idx].depth = depth() + 1;

    if(!wasLeaf)
        ++(tree->numLeafs);
    if(child.depth() > tree->depth)
        tree->depth = child.depth();

    if(n && n->idx >= i) {
        Debug{} << "Adjusting node";
        ++(n->idx);
    }
    return child;
}

Node Node::addRightChild(Node* n) { return addChild(false, n); }

Node Node::addLeftChild(Node* n) { return addChild(true, n); }

void Node::initializePhasefieldFromParent() {
    Mesh* mesh = tree->mesh;
    mesh->requireFaceInformation();

    // @TODO weights only need to be computed for the a single node..
    tree->computeWeightsOfAncestorsOfLevel(depth());
    auto weights = temporary();
    auto phase = phasefield();

    FaceData<char> thresholded{NoInit, mesh->faceCount()};
    double totalArea = 0.;
    Face start{Invalid, mesh};
    for(Face f : mesh->faces()) {
        double v = 0;
        for(Vertex vertex : f.vertices())
            v += weights[vertex];
        if(v/3. > 1e-2) {
            totalArea += f.area();
            thresholded[f] = 1;
            start = f;
        } else
            thresholded[f] = 0;
    }

    if(!!start) {
        Bfs<Face> bfs{*mesh};
        bfs.setSource(start);
        double area = 0;
        Face f;
        while(bfs.step(f)) {
            area += f.area();
            if(area > totalArea*0.5)
                break;
        }

        for(Face face : mesh->faces()) {
            if(thresholded[face] > 0.) {
                for(Vertex v : face.vertices()) {
                    phase[v] = bfs.visited(face) ? 1. : -1.;
                }
            }
        }
    } else Debug{} << "Could not find any face for which weights are positive";

}

void Node::splitAndInitialize(Node* n) {
    Node right;
    if(!hasRightChild()) {
        right = addRightChild(n);
    } else right = rightChild();
    right.initializePhasefieldFromParent();

    Node left;
    if(!hasLeftChild()) {
        left = addLeftChild(n);
    } else left = leftChild();
    left.initializePhasefieldFromParent();
}

double Node::integrateWeight(Mesh& mesh) const {
    mesh.requireIntegralOperator();
    auto weights = temporary();
    double sum = 0;
    for(Vertex v : mesh.vertices())
        sum += mesh.integral[v]*weights[v];

    return sum;
}


Debug& operator<<(Debug& debug, Node const& n) {
    std::string leftChild = n.leftChild().isValid() ? std::to_string(n.leftChild().idx) : "Invalid";
    std::string rightChild = n.rightChild().isValid() ? std::to_string(n.rightChild().idx) : "Invalid";
    std::string parent = n.parent().isValid() ? std::to_string(n.parent().idx) : "Invalid";
    std::string formatted = Cr::Utility::formatString(
            "Node (idx = {}, depth = {}, left child = {}, right child = {}, parent  = {})", n.idx, n.depth(), leftChild,
            rightChild, parent);
    debug << formatted.c_str();
    return debug;
}

Tree::Tree(Mesh& m) : mesh(&m), nodeData(1), numLeafs(1) { mesh->requireFaceInformation(); }

void Tree::update() {
    size_t n = mesh->vertexCount();
    size_t min = Math::min(n, phasefieldData.size()/nodeCount());

    Array<double> data(DirectInit, n*nodeCount(), 0.);

    if(min) {
        StridedArrayView2D<double> view{data, {nodeData.size(), n}};
        size_t i = 0;
        for(Node node : nodes()) {
            Cr::Utility::copy(node.phasefield().prefix(min), view[i++].prefix(min));
        }
    }

    arrayResize(tempsData, data.size());
    phasefieldData = std::move(data);

}

Node Tree::insertNodeAtIndex(size_t idx) {
    size_t n = nodeCount();
    size_t m = vertexCount();

    /* reset pointers*/
    for(NodeData& nd : nodeData) {
        if(nd.parent != Invalid && nd.parent >= idx) ++nd.parent;
        if(nd.leftChild != Invalid && nd.leftChild >= idx) ++nd.leftChild;
        if(nd.rightChild != Invalid && nd.rightChild >= idx) ++nd.rightChild;
    }

    arrayResize(tempsData, NoInit, m*(n + 1));
    arrayResize(phasefieldData, NoInit, m*(n + 1));
    arrayResize(nodeData, NoInit, n + 1);

    double* ptrTemp = tempsData.data();
    double* ptrPhase = phasefieldData.data();
    NodeData* nodeDataPtr = nodeData.data();

    /* move everything back by one */
    memmove(ptrTemp + m*(idx + 1), ptrTemp + m*idx, (n - idx)*m*sizeof(double));
    memmove(ptrPhase + m*(idx + 1), ptrPhase + m*idx, (n - idx)*m*sizeof(double));
    memmove(nodeDataPtr + idx + 1, nodeDataPtr + idx, (n - idx)*sizeof(NodeData));

    nodeData[idx].depth = Invalid;
    nodeData[idx].leftChild = Invalid;
    nodeData[idx].rightChild = Invalid;
    nodeData[idx].parent = Invalid;

    Node node{idx, this};
    for(double& x : node.phasefield()) x = 0.;

    return node;
}


ArrayView<double> Tree::level(size_t d) {
    size_t begin = levelStartIndex(d);
    size_t end = levelStartIndex(d + 1);
    return phasefieldData.slice(begin*vertexCount(), end*vertexCount());
}

void Tree::serialize(Array<char>& data) const {
    size_t dataSize = phasefieldData.size();
    arrayAppend(data, {(char*) &dataSize, sizeof(size_t)});
    arrayAppend(data, arrayCast<char>(phasefieldData));
    size_t nodesSize = nodeData.size();
    arrayAppend(data, {(char*) &nodesSize, sizeof(size_t)});
    arrayAppend(data, arrayCast<char>(nodeData));
}

Tree Tree::deserialize(ArrayView<const char> const& data, Mesh& m) {
    Tree t{m};
    char const* pc = data;
    auto dataSize = deserializeTrivial<size_t>(pc);
    arrayResize(t.phasefieldData, NoInit, dataSize);
    arrayResize(t.tempsData, NoInit, dataSize);
    memcpy(t.phasefieldData.data(), pc, sizeof(double)*dataSize);
    pc += sizeof(double)*dataSize;

    auto nodesSize = deserializeTrivial<size_t>(pc);
    arrayResize(t.nodeData, NoInit, nodesSize);
    memcpy(t.nodeData.data(), pc, sizeof(NodeData)*nodesSize);

    pc += sizeof(Node)*nodesSize;

    t.numLeafs = 0;
    t.depth = 0;
    for(Node node : t.nodes()) {
        if(node.isLeaf()) t.numLeafs++;
        t.depth = Math::max(t.depth, node.depth());
    }

    return t;
}

size_t Tree::levelStartIndex(size_t level) {
    for(size_t i = 0; i < nodeData.size(); ++i) {
        if(nodeData[i].depth == level) return i;
    }
    return nodeData.size();
}

Range<LeafIterator> Tree::leafs() {
    LeafIterator b{0, this};
    if(!b.node.isLeaf()) ++b;
    return {b, {nodeCount(), this}};
}

Range<HorizontalNodeIterator> Tree::nodesOnLevel(size_t l) {
    size_t b = levelStartIndex(l);
    size_t e = levelStartIndex(l + 1);
    return {{b, this},
            {e, this}};
}

Range<HorizontalNodeIterator> Tree::nodesBelowLevel(size_t l) { return {{0,                      this},
                                                                        {levelStartIndex(l + 1), this}};
}

Range<HorizontalNodeIterator> Tree::nodes() { return {{0,           this},
                                                      {nodeCount(), this}};
}

Range<InternalNodeIterator> Tree::internalNodes() {
    InternalNodeIterator b{0, this};
    if(b.node.isLeaf()) ++b;
    return {b, {nodeCount(), this}};
}

Array<Node> Tree::ancestorsOfLevel(size_t l) {
    Array<Node> ancestors;
    for(Node node : nodesOnLevel(l)) {
        while(node.parent()) { /* walk up the tree */
            arrayAppend(ancestors, node.parent());
            node = node.parent();
        }
    }

    std::sort(ancestors.begin(), ancestors.end());
    Node* it = std::unique(ancestors.begin(), ancestors.end());
    arrayResize(ancestors, it - ancestors.begin());
    return ancestors;
}

void Tree::computeWeightsOfAncestorsOfLevel(size_t l) {
    size_t n = vertexCount();
    SmootherStep smoothStep;

    for(double& w : root().temporary()) w = 1.;

    for(Node node : ancestorsOfLevel(l)) {
        auto weights = node.temporary();
        for(size_t i = 0; i < n; ++i) {
            if(node.hasLeftChild()) {
                Node leftChild = node.leftChild();
                leftChild.temporary()[i] = smoothStep.eval(node.phasefield()[i])*weights[i];
            }
            if(node.hasRightChild()) {
                Node rightChild = node.rightChild();
                rightChild.temporary()[i] = smoothStep.eval(-node.phasefield()[i])*weights[i];
            }
        }
    }
}

void Tree::reset() {
    size_t n = vertexCount();
    arrayResize(nodeData, 1);
    arrayResize(phasefieldData, n);
    arrayResize(tempsData, n);
    depth = 0;
    numLeafs = 1;
    nodeData[0].rightChild = Invalid;
    nodeData[0].leftChild = Invalid;
}

void Tree::computeLeafWeights() {
    size_t n = vertexCount();
    SmootherStep smoothStep;

    for(double& w : root().temporary()) w = 1.;

    for(Node node : nodes()) {
        if(node.isLeaf()) continue;

        auto weights = node.temporary();
        for(size_t i = 0; i < n; ++i) {
            if(node.hasLeftChild()) {
                Node leftChild = node.leftChild();
                leftChild.temporary()[i] = smoothStep.eval(node.phasefield()[i])*weights[i];
            }
            if(node.hasRightChild()) {
                Node rightChild = node.rightChild();
                rightChild.temporary()[i] = smoothStep.eval(-node.phasefield()[i])*weights[i];
            }
        }
    }
}

StridedArrayView2D<double> Tree::phasefields() { return {phasefieldData, {nodeCount(), vertexCount()}}; }

}