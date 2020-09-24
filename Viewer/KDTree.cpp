//
// Created by janos on 9/21/20.
//

#include "KDTree.h"
#include "StlAlgorithm.h"
#include "Mesh.h"

#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/StridedArrayView.h>

#include <Magnum/Math/FunctionsBatch.h>
#include <Magnum/Math/Functions.h>

#include <string>

namespace Phasefield {

KDTree::KDTree(Mesh& mesh) : MeshFeature(mesh, false), m_mesh(&mesh) { updateImpl(); }

void KDTree::update() { updateImpl(); }

Vector3 const& KDTree::point(NodeHandle handle) const { return points()[m_nodes[handle].pointIndex]; }

size_t KDTree::leftChild(NodeHandle handle) const { return m_nodes[handle].leftChild; }

size_t KDTree::rightChild(NodeHandle handle) const { return m_nodes[handle].rightChild; }

size_t KDTree::root() const { return m_root; }

size_t KDTree::constructTree(size_t* begin, size_t* end, size_t depth, size_t& size) {
    if(end <= begin) return Invalid;

    auto cd = depth%3;
    auto n = begin + (end - begin)/2;

    std::nth_element(begin, n, end,
                     [&](size_t id1, size_t id2) { return points()[id1][cd] < points()[id2][cd]; });

    auto handle = size++;
    auto& node = m_nodes[handle];
    node.pointIndex = *n;
    node.leftChild = constructTree(begin, n, depth + 1, size);
    node.rightChild = constructTree(n + 1, end, depth + 1, size);
    return handle;
}

void KDTree::recurse(Vector3 const& q, size_t nodeId, size_t cd, NNResult& result) {
    if(nodeId == Invalid) return;
    auto const& node = m_nodes[nodeId];
    auto const& p = points()[node.pointIndex];
    /* prune by computing distance to bounding box, @todo does this actually speed things up */
    //if(bbDistanceSq(q, bb) > result.distanceSquared) return;
    auto distSq = (p - q).dot();
    if(distSq < result.distanceSquared) {
        result.distanceSquared = distSq;
        result.pointIndex = node.pointIndex;
    }
    auto nextCd = (cd + 1)%3;
    if(q[cd] < p[cd]) { /* q is closer to left child */
        recurse(q, node.leftChild, nextCd, result);
        /* prune by computing distance to splitting plane */
        if(p[cd] - q[cd] < result.distanceSquared)
            recurse(q, node.rightChild, nextCd, result);
    } else { /* q is closer to right child */
        recurse(q, node.rightChild, nextCd, result);
        if(q[cd] - p[cd] < result.distanceSquared)
            recurse(q, node.leftChild, nextCd, result);
    }
}

KDTree::NNResult KDTree::nearestNeighbor(Vector3 const& queryPoint) {
    NNResult result;
    recurse(queryPoint, m_root, 0, result);
    return result;
}

void KDTree::updateImpl() {
    arrayResize(m_nodes, NoInit, vertexCount());
    Array<size_t> nodeIds(NoInit, vertexCount());
    for(size_t i = 0; i < nodeIds.size(); ++i)
        nodeIds[i] = i;
    size_t size = 0;
    m_root = constructTree(nodeIds.begin(), nodeIds.end(), 0, size);
}

size_t KDTree::vertexCount() const {
    return m_mesh->vertexCount();
}

StridedArrayView1D<const Vector3> KDTree::points() const { return m_mesh->positions(); }

void formatTree(
        Magnum::Debug& debug,
        std::string const& prefix,
        KDTree const& tree,
        KDTree::NodeHandle handle,
        bool isLeft,
        const char* arrow) {
    if(handle != Invalid) {
        debug << prefix.c_str() << arrow << tree.point(handle) << "\n";
        // enter the next tree level - left and right branch

        const char* leftChildArrow = "├──";
        if(tree.rightChild(handle) == Invalid)
            leftChildArrow = "└──";
        formatTree(debug, prefix + (isLeft ? " │   " : "     "), tree, tree.leftChild(handle), true, leftChildArrow);
        formatTree(debug, prefix + (isLeft ? " │   " : "     "), tree, tree.rightChild(handle), false, "└──");
    }
}

Magnum::Debug& operator<<(Magnum::Debug& debug, KDTree const& tree) {
    formatTree(debug, "", tree, tree.root(), false, "└──");
    return debug;
}

}