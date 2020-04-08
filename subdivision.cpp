//
// Created by janos on 03.04.20.
//

#include "subdivision.hpp"
#include "upload.hpp"

#include <Magnum/MeshTools/Subdivide.h>
#include <Magnum/MeshTools/Interleave.h>

#include <imgui.h>

using namespace Magnum;
using namespace Corrade;

void Subdivision::subdivide(int numSubdivisions) {

    Containers::Array<UnsignedInt> faces;
    Containers::Array<Vector3> vertices;
    {
        std::lock_guard l(m_phasefieldData.mutex);
        faces = m_phasefieldData.original.indicesAsArray();
        vertices = m_phasefieldData.original.positions3DAsArray();
    }

    while(numSubdivisions--)
        MeshTools::subdivide(faces,vertices,[](auto const& v1, auto const& v2){ return .5 * (v1 + v2); });

    Trade::MeshAttributeData vertexData{Trade::MeshAttribute::Position, VertexFormat::Vector3, vertices};

    Trade::MeshData meshData{MeshPrimitive::Triangles,
                         Trade::DataFlag::Mutable, faces, Trade::MeshIndexData{faces},
                         Trade::DataFlag::Mutable, vertices, {vertexData}};

    {
        std::lock_guard l(m_phasefieldData.mutex);
        auto& pd = m_phasefieldData;
        pd.meshData = preprocess(meshData, CompileFlag::GenerateSmoothNormals|CompileFlag::AddTextureCoordinates);
        pd.status = PhasefieldData::Status::NewMesh;

        pd.V = std::move(vertices);
        pd.F = std::move(faces);
    }
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
