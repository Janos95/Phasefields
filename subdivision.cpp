//
// Created by janos on 03.04.20.
//

#include "subdivision.hpp"
#include "upload.hpp"

#include <Corrade/Utility/Algorithms.h>
#include <Magnum/MeshTools/Subdivide.h>
#include <Magnum/MeshTools/RemoveDuplicates.h>
#include <Magnum/MeshTools/Interleave.h>

#include <imgui.h>
#include <fmt/core.h>

using namespace Magnum;
using namespace Corrade;

void Subdivision::subdivide(int numSubdivisions) {
    auto& pd = m_phasefieldData;

    Trade::MeshData* data;
    Containers::ArrayView<Float> phasefield;
    if(numSubdivisions >= m_lastNumSubdivision){
        numSubdivisions -= m_lastNumSubdivision;
        data = &pd.meshData;
        phasefield = pd.phasefield;
    } else {
        data = &pd.original;
        phasefield = Containers::ArrayView<Float>{pd.phasefield, data->vertexCount()};
    }

    auto positions = data->attribute<Vector3>(Trade::MeshAttribute::Position);
    auto indices = data->indicesAsArray();
    auto vertexData = MeshTools::interleave(positions, phasefield);

    auto indicesSizeCurrent = indices.size();
    auto verticesSizeCurrent = positions.size();
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
    fmt::print("Old size was {}, new size after removing is {}\n", verticesSize, sizeAfterRemoving);

    Trade::MeshAttributeData positionData{Trade::MeshAttribute::Position, VertexFormat::Vector3, {vertexData, &vertices[0].xyz(), sizeAfterRemoving, sizeof(Vector4)}};
    Trade::MeshData meshData{MeshPrimitive::Triangles,
                         Trade::DataFlag{}, indices, Trade::MeshIndexData{indices},
                         Trade::DataFlag{}, vertices, {positionData}};

    pd.meshData = preprocess(meshData, CompileFlag::GenerateSmoothNormals|CompileFlag::AddTextureCoordinates);
    pd.V = pd.meshData.positions3DAsArray();
    pd.F = pd.meshData.indicesAsArray();

    Containers::arrayResize(pd.phasefield, sizeAfterRemoving);
    auto textureView = pd.meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);
    for (int j = 0; j < sizeAfterRemoving; ++j) {
        auto u = vertices[j].w();
        pd.phasefield[j] = u;
        textureView[j].x() = .5f * (u + 1.f);
    }

    pd.status = PhasefieldData::Status::Subdivided;
}

void Subdivision::drawImGui(){
    if (ImGui::TreeNode("Subdivisions"))
    {
        ImGui::Text("Currently we have %d vertices and %d faces", m_numVertices, m_numFaces);
        static int numSubdivisions = 1;
        constexpr int step = 1;
        ImGui::InputScalar("Number of Sibdivision (wrt. orignal mesh)", ImGuiDataType_U32, &numSubdivisions, &step, nullptr, "%d");
        if(ImGui::Button("do subdivision")){
            subdivide(numSubdivisions);
            m_lastNumSubdivision= numSubdivisions;
            m_numVertices = m_phasefieldData.meshData.vertexCount();
            m_numFaces = m_phasefieldData.meshData.indexCount() / 3;
        }
        ImGui::TreePop();
    }
}
