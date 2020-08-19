//
// Created by janos on 16.06.20.
//

#ifndef KDTREE_KDTREE_H
#define KDTREE_KDTREE_H

#include <Magnum/Magnum.h>
#include <Magnum/Math/Vector3.h>
#include <Magnum/Math/Distance.h>
#include <Magnum/Math/Range.h>
#include <Magnum/Math/Functions.h>
#include <Magnum/Math/FunctionsBatch.h>
#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/StridedArrayView.h>

#include <algorithm>
#include <string>

using namespace Magnum;
using namespace Corrade;

enum class Side : UnsignedByte {
    Left,
    Right
};

template<Side side, UnsignedInt Size, class T>
Math::Range<Size, T> trim(
        Math::Range<Size, T> range,
        UnsignedInt dim,
        typename Math::Range<Size,T>::VectorType const& p)
{
    if constexpr(side == Side::Right) range.min()[dim] = p[dim];
    else range.max()[dim] = p[dim];
    return range;
}

template<class Vector>
typename Vector::Type bbDistanceSq(Vector const& p, Math::Range<Vector::Size, typename Vector::Type> const& range){
    using Type = typename Vector::Type;
    Type dsq{0};
    for (int j = 0; j < Vector::Size; ++j){
        auto lower = range.min()[j];
        auto upper = range.max()[j];
        if(p[j] < lower)
            dsq += Math::pow(lower - p[j], Type{2.});
        else if(p[j] > upper)
            dsq += Math::pow(p[j] - upper, Type{2.});
    }
    return dsq;
}

template<class Vector>
class KDTree {
public:

    struct Node {
        constexpr static int32_t Invalid = -1;
        int32_t leftChild;
        int32_t rightChild;
        int32_t pointIndex;
    };

    enum {
        Size = Vector::Size
    };
    using Type = typename Vector::Type;
    using Range = Math::Range<Size, Type>;
    using NodeHandle = int32_t;

    KDTree() = default;

    explicit KDTree(Containers::StridedArrayView1D<const Vector> const& points) :
            m_nodes(Containers::NoInit, points.size()),
            m_points(points),
            m_bb(Math::minmax(points))
    {
        Containers::Array<UnsignedInt> nodeIds(Containers::NoInit, points.size());
        for (UnsignedInt i = 0; i < nodeIds.size(); ++i)
            nodeIds[i] = i;
        uint32_t size = 0;
        m_root = constructTree(nodeIds.begin(), nodeIds.end(), 0, size);
    }

    auto& point(NodeHandle handle) const { return m_points[m_nodes[handle].pointIndex]; }
    auto leftChild(NodeHandle handle) const { return m_nodes[handle].leftChild; }
    auto rightChild(NodeHandle handle) const { return m_nodes[handle].rightChild; }
    auto root() const { return m_root; }
    auto bb() const { return m_bb; }

    struct NNResult{
        int pointIndex;
        Type distanceSquared = Math::Constants<Type>::inf();
    };

    /**
     * @param queryPoint the point for which to perform
     * the nearest neighbor query.
     * @return struct containing index corresponding to nearst neighbor and
     * squared distance
     */
    NNResult nearestNeighbor(Vector const& queryPoint){
        NNResult result;
        recurse(queryPoint, m_root, 0, m_bb, result);
        return result;
    }

private:

    uint32_t constructTree(uint32_t* begin, uint32_t* end, uint32_t depth, uint32_t& size){
        if(end <= begin) return Node::Invalid;

        auto cd = depth % Size;
        auto n = begin + (end - begin)/2;

        std::nth_element(begin, n, end,
                         [&](uint32_t id1, uint32_t id2){ return m_points[id1][cd] < m_points[id2][cd]; });

        auto handle = size++;
        auto& node = m_nodes[handle];
        node.pointIndex = *n;
        node.leftChild = constructTree(begin, n, depth + 1, size);
        node.rightChild = constructTree(n + 1, end, depth + 1, size);
        return handle;
    }

    void recurse(Vector const& q, int nodeId, int cd, Range const& bb, NNResult& result){
        if(nodeId == Node::Invalid) return;
        auto const& node = m_nodes[nodeId];
        auto const& p = m_points[node.pointIndex];
        /* prune by computing distance to bounding box, @todo does this actually speed things up */
        if(bbDistanceSq(q, bb) > result.distanceSquared) return;
        auto distSq = (p - q).dot();
        if(distSq < result.distanceSquared){
            result.distanceSquared = distSq;
            result.pointIndex = node.pointIndex;
        }
        auto nextCd = (cd+1) % Size;
        if(q[cd] < p[cd]){ /* q is closer to left child */
            recurse(q, node.leftChild, nextCd, trim<Side::Left>(bb, cd, p), result);
            /* prune by computing distance to splitting plane */
            if(p[cd] - q[cd] < result.distanceSquared)
                recurse(q, node.rightChild, nextCd, trim<Side::Right>(bb, cd, p), result);
        } else { /* q is closer to right child */
            recurse(q, node.rightChild, nextCd, trim<Side::Right>(bb, cd, p), result);
            if(q[cd] - p[cd] < result.distanceSquared)
                recurse(q, node.leftChild, nextCd, trim<Side::Left>(bb, cd, p), result);
        }
    }

    int m_root;
    Range m_bb;
    Containers::Array<Node> m_nodes;
    Containers::StridedArrayView1D<const Vector> m_points;
};

/* deduction guide */
template<class Vector>
KDTree(Containers::Array<Vector> const& points) -> KDTree<Vector>;

template<class T>
void formatTree(
        Magnum::Debug& debug,
        std::string const& prefix,
        KDTree<T> const& tree,
        typename KDTree<T>::NodeHandle handle,
        bool isLeft,
        const char* arrow)
{
    if(handle != KDTree<T>::Node::Invalid){
        debug << prefix.c_str() << arrow << tree.point(handle) << "\n";
        // enter the next tree level - left and right branch

        const char* leftChildArrow = "├──";
        if(tree.rightChild(handle) == KDTree<T>::Node::Invalid)
            leftChildArrow = "└──";
        formatTree(debug, prefix + (isLeft ? " │   " : "     "), tree, tree.leftChild(handle), true, leftChildArrow);
        formatTree(debug, prefix + (isLeft ? " │   " : "     "), tree, tree.rightChild(handle), false, "└──");
    }
}

template<class T>
Magnum::Debug& operator<<(Magnum::Debug& debug, KDTree<T> const& tree){
    formatTree(debug, "", tree, tree.root(), false, "└──");
    return debug;
}

#endif //KDTREE_KDTREE_H