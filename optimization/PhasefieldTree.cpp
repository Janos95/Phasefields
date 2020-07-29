//
// Created by janos on 7/27/20.
//

#include "PhasefieldTree.h"

#include <Corrade/Utility/Algorithms.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Utility/MurmurHash2.h>

#include <Magnum/Math/Vector4.h>
#include <Magnum/Magnum.h>

#include <cstdio>
#include <unordered_map>

using namespace Magnum;
using namespace Corrade;

namespace {

struct Edge {
    Edge(Mg::UnsignedInt v1_, Mg::UnsignedInt v2_) : v1(std::min(v1_, v2_)), v2(std::max(v1_, v2_)) {}

    Mg::UnsignedInt v1, v2;

    auto operator<=>(const Edge&) const = default;
};

struct EdgeHash {
    std::size_t operator()(Edge const& e) const {
        return *((std::size_t*) (Utility::MurmurHash2{}((char*) &e, sizeof(UnsignedLong)).byteArray()));
    }
};

}

void PhasefieldTree::subdivide(Containers::Array<Mg::UnsignedInt>& indices, Containers::Array<Mg::Vector3d>& vertices){

    std::unordered_map<Edge, int, EdgeHash> edgeMap(indices.size());

    const auto pfs = phasefields();
    Containers::Array<Containers::Array<double>> interpolated(phasefieldCount());

    const auto vertexOffset = vertices.size();
    int edgeCount = 0;

    const std::size_t indexCount = indices.size();

    for(std::size_t i = 0; i < indexCount; i += 3) {
        UnsignedInt newVertices[3];
        for(int j = 0; j < 3; ++j) {
            const auto v1 = indices[i + j], v2 = indices[i + (j + 1)%3];
            const auto[it, inserted] = edgeMap.try_emplace(Edge{v1, v2}, vertexOffset + edgeCount);
            newVertices[j] = it->second;
            if(inserted) {
                Containers::arrayAppend(vertices, 0.5*(vertices[v1] + vertices[v2]));
                for(std::size_t k = 0; k < phasefieldCount(); ++k) {
                    Containers::arrayAppend(interpolated[k], 0.5*(pfs[j][v1] + pfs[j][v2]));
                }
                ++edgeCount;
            }
        }

        /*
            Add three new faces (0, 1, 3) and update original (2)

                          orig 0
                          /   \
                         /  0  \
                        /       \
                    new 0 ----- new 2
                    /   \       /  \
                   /  1  \  2  / 3  \
                  /       \   /      \
             orig 1 ----- new 1 ---- orig 2
        */
        Containers::arrayAppend(indices, {
                indices[i], newVertices[0], newVertices[2],
                newVertices[0], indices[i + 1], newVertices[1],
                newVertices[2], newVertices[1], indices[i + 2]});

        for(std::size_t j = 0; j != 3; ++j)
            indices[i + j] = newVertices[j];
    }

    /* flatten everything */
    Containers::Array<double> data(phasefieldCount()*vertices.size());
    Containers::StridedArrayView2D<double> phasefieldsNew(phasefieldData, {phasefieldCount(), vertices.size()});
    for(std::size_t i = 0; i < phasefieldCount(); ++i) {
        Utility::copy(pfs[i], phasefieldsNew[i].prefix(phasefieldSize));
        Utility::copy(interpolated[i], phasefieldsNew[i].slice(phasefieldSize, vertices.size()));
    }

    /* mutate tree structure */
    phasefieldData = std::move(data);
    Containers::arrayResize(tempsData, phasefieldData.size());
    phasefieldSize = vertices.size();
}


