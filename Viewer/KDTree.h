//
// Created by janos on 16.06.20.
//

#pragma once

#include "Types.h"
#include "Surface.h"
#include "MeshFeature.h"

#include <Corrade/Containers/Array.h>

#include <Magnum/Math/Range.h>
#include <Magnum/Math/Constants.h>

namespace Phasefield {

namespace Mg = Magnum;

enum class Side : Mg::UnsignedByte {
    Left,
    Right
};

class KDTree : MeshFeature {
public:

    struct Node {
        size_t leftChild;
        size_t rightChild;
        size_t pointIndex;
    };

    using NodeHandle = size_t;

    KDTree() = default;

    explicit KDTree(Mesh& mesh);

    [[nodiscard]] Vector3 const& point(NodeHandle handle) const;

    [[nodiscard]] size_t leftChild(NodeHandle handle) const;

    [[nodiscard]] size_t rightChild(NodeHandle handle) const;

    [[nodiscard]] size_t root() const;

    void updateImpl();
    void update() override;

    FEATURE_NAME("KDTree")

    struct NNResult {
        size_t pointIndex;
        float distanceSquared = Mg::Math::Constants<float>::inf();
    };

    /**
     * @param queryPoint the point for which to perform
     * the nearest neighbor query.
     * @return struct containing index corresponding to nearst neighbor and
     * squared distance
     */
    NNResult nearestNeighbor(Vector3 const& queryPoint);

private:

    [[nodiscard]] StridedArrayView1D<const Vector3> points() const;

    [[nodiscard]] size_t vertexCount() const;

    size_t constructTree(size_t* begin, size_t* end, size_t depth, size_t& size);

    void recurse(Vector3 const& q, size_t nodeId, size_t cd, NNResult& result);

    Mesh const* m_mesh = nullptr;
    size_t m_root;
    Array<Node> m_nodes;
};

Magnum::Debug& operator<<(Magnum::Debug& debug, KDTree const& tree);

}
