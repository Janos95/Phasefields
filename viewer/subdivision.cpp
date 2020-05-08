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


using namespace Magnum;
using namespace Corrade;

void subdivide(
        std::uint32_t numSubdivisions,
        Trade::MeshData const& original,
        Containers::Array<Double>& phasefield,
        Trade::MeshData& meshData) {

    Containers::Array<Float> phasefieldF(Containers::NoInit, original.vertexCount());
    for(auto i = 0; i < original.vertexCount(); ++i) phasefieldF[i] = phasefield[i];
    auto indices = original.indicesAsArray();
    auto vertexData = MeshTools::interleave(original.attribute<Vector3>(Trade::MeshAttribute::Position), phasefieldF);

    auto indicesSizeCurrent = indices.size();
    auto verticesSizeCurrent = original.vertexCount();
    auto indicesSize = indicesSizeCurrent * std::pow(4, numSubdivisions);
    auto verticesSize = verticesSizeCurrent + (indicesSize - indicesSizeCurrent)/3;

    Containers::arrayResize(indices,indicesSize);
    Containers::arrayResize(vertexData,verticesSize*sizeof(Vector4));

    auto vertices = Containers::arrayCast<Vector4>(vertexData);

    while(numSubdivisions--){
        Containers::ArrayView<UnsignedInt> indicesView{indices.data(), indicesSizeCurrent * 4};
        Containers::StridedArrayView1D<Vector4> verticesView{vertices, verticesSizeCurrent + indicesSizeCurrent};

        MeshTools::subdivideInPlace(indicesView, verticesView, [](auto& v1, auto& v2){return .5f * (v1 + v2); });

        verticesSizeCurrent += indicesSizeCurrent;
        indicesSizeCurrent *= 4;
    }

    auto sizeAfterRemoving = MeshTools::removeDuplicatesIndexedInPlace<UnsignedInt, Vector4>(indices, vertices);
    printf("Old size was %d, new size after removing is %d\n", (int)verticesSize, (int)sizeAfterRemoving);

    Containers::StridedArrayView1D<Vector3> stridedView{vertexData, &vertices[0].xyz(), sizeAfterRemoving, sizeof(Vector4)};
    Trade::MeshAttributeData positionData{Trade::MeshAttribute::Position, stridedView};
    Trade::MeshData md{MeshPrimitive::Triangles,
                         {}, indices, Trade::MeshIndexData{indices},
                         {}, vertexData, {positionData}};

    //update everything
    meshData = preprocess(md, CompileFlag::GenerateSmoothNormals|CompileFlag::AddTextureCoordinates);
    Containers::arrayResize(phasefield, sizeAfterRemoving);
    for (int j = 0; j < sizeAfterRemoving; ++j) {
        auto u = vertices[j].w();
        phasefield[j] = u;
    }
}

