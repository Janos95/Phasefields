//
// Created by janos on 03.04.20.
//

#include "subdivision.hpp"
#include "upload.hpp"

#include <Corrade/Utility/Algorithms.h>
#include <Magnum/MeshTools/Subdivide.h>
#include <Magnum/MeshTools/RemoveDuplicates.h>

#include <imgui.h>
#include <fmt/core.h>

using namespace Magnum;
using namespace Corrade;

void Subdivision::subdivide(int numSubdivisions) {
    auto faces = m_phasefieldData.original.indicesAsArray();
    auto vertices = m_phasefieldData.original.positions3DAsArray();

    /**
     * @todo this can be done with only one (re-)allocation, by resizing upfront.
     */
    while(numSubdivisions--){
        auto linearInterpolation = [](auto const& v1, auto const& v2){ return .5 * (v1 + v2); };
        MeshTools::subdivide(faces,vertices,linearInterpolation);
    }

    auto oldSize = vertices.size();
    auto newSize = MeshTools::removeDuplicatesIndexedInPlace<UnsignedInt, Vector3>(faces, vertices);
    fmt::print("Old size was {}, new size after removing is {}\n", oldSize, newSize);

    Trade::MeshAttributeData vertexData{Trade::MeshAttribute::Position, VertexFormat::Vector3, vertices};

    Trade::MeshData meshData{MeshPrimitive::Triangles,
                         Trade::DataFlag::Mutable, faces, Trade::MeshIndexData{faces},
                         Trade::DataFlag::Mutable, vertices, {vertexData}};

    auto& pd = m_phasefieldData;
    pd.meshData = preprocess(meshData, CompileFlag::GenerateSmoothNormals|CompileFlag::AddTextureCoordinates);
    pd.status = PhasefieldData::Status::NewMesh;

    pd.V = std::move(vertices);
    pd.F = std::move(faces);
    /**
     * @todo interpolate the phasefield instead of setting it to zero.
     * But this should handle the case when the number of subdivision is less
     * than what was used before
     */
    Containers::arrayResize(pd.phasefield, pd.V.size());
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
        }
        ImGui::TreePop();
    }
}
