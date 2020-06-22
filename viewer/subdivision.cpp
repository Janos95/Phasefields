//
// Created by janos on 03.04.20.
//

#include "subdivision.hpp"
#include "upload.hpp"

#include <Corrade/Utility/Algorithms.h>

#include <Magnum/MeshTools/Subdivide.h>
#include <Magnum/MeshTools/RemoveDuplicates.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/Math/Vector4.h>
#include <Magnum/Magnum.h>


using namespace Magnum;
using namespace Corrade;

void subdivide(
        std::uint32_t numSubdivisions,
        Containers::Array<UnsignedInt>& indices,
        Containers::Array<Vector3d>& vertices,
        Containers::Array<Double>& phasefield) {
    auto vertexData = MeshTools::interleave(vertices, phasefield);

    auto indicesSizeCurrent = indices.size();
    auto verticesSizeCurrent = vertices.size();
    auto indicesSize = indicesSizeCurrent * std::pow(4, numSubdivisions);
    auto verticesSize = verticesSizeCurrent + (indicesSize - indicesSizeCurrent)/3;

    Containers::arrayResize(indices,indicesSize);
    Containers::arrayResize(vertexData,verticesSize*sizeof(Vector4d));

    auto dataView = Containers::arrayCast<Vector4d>(vertexData);

    while(numSubdivisions--){
        Containers::ArrayView<UnsignedInt> indicesView{indices.data(), indicesSizeCurrent * 4};
        Containers::StridedArrayView1D<Vector4d> verticesView{dataView, verticesSizeCurrent + indicesSizeCurrent};
        MeshTools::subdivideInPlace(indicesView, verticesView, [](auto& v1, auto& v2){return .5f * (v1 + v2); });

        verticesSizeCurrent += indicesSizeCurrent;
        indicesSizeCurrent *= 4;
    }

    auto sizeAfterRemoving = MeshTools::removeDuplicatesFuzzyIndexedInPlace(indices, Containers::arrayCast<2, double>(dataView));
    printf("Old size was %d, new size after removing is %d\n", (int)verticesSize, (int)sizeAfterRemoving);

    Containers::arrayResize(phasefield, sizeAfterRemoving);
    Containers::arrayResize(vertices, sizeAfterRemoving);
    for (int j = 0; j < sizeAfterRemoving; ++j) {
        phasefield[j] = dataView[j].w();
        vertices[j] = dataView[j].xyz();
    }
}

