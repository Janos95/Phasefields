//
// Created by janos on 04.04.20.
//

#include "primitives.hpp"

#include <Corrade/Utility/Algorithms.h>

#include <Magnum/Primitives/Capsule.h>
#include <Magnum/MeshTools/Interleave.h>

#include <folly/sorted_vector_types.h>

#include <imgui.h>

using ComboElement = LoadPrimitives::ComboElement;

namespace {
    struct CompareByName{
        bool operator()(ComboElement const& e1, ComboElement const& e2){
            return e1.name < e2.name;
        }
    };
}

auto makeComboMap(){
   std::vector<ComboElement> map;
   map.push_back(ComboElement{"U shaped square", PrimitiveType::U, std::make_unique<LoadPrimitives::UOptions>()});
   map.push_back(ComboElement{"Capsule", PrimitiveType::Capsule, std::make_unique<LoadPrimitives::CapsuleOptions>()});
   std::sort(map.begin(), map.end(), CompareByName{});
   return map;
}


void LoadPrimitives::drawImGui(){
    if (ImGui::TreeNode("Primitives"))
    {
        static auto map = makeComboMap();
        auto current = map.end();

        if (ImGui::BeginCombo("##combo", current != map.end() ? current->name.c_str() : nullptr)) {
            for (auto it = map.begin(); it < map.end(); ++it){
                bool isSelected = (current == it); // You can store your selection however you want, outside or inside your objects
                if (ImGui::Selectable(it->name.c_str(), isSelected))
                    current = it;
                if (isSelected)
                    ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
            }
            //@todo imgui display for options
            ImGui::EndCombo();
        }

        if(ImGui::Button("Load") && current != map.end()){
            if(!m_phasefieldData)
            load(*current);
        }

        ImGui::TreePop();
    }
}


Trade::MeshData uShapedSquare(float width, float height, float innerWidth, float innerHeight){
    float a = .5f * (width - innerWidth);
    float b = a + innerWidth;

    Containers::Array<Vector3> vertices{
        Containers::InPlaceInit,
        {
            {.0f, .0f, .0f},
            {a, .0f, .0f},
            {a, innerHeight, .0f},
            {.0f, innerHeight, .0f},
            {.0f, height, .0f},
            {a, height, .0f},
            {b, height, .0f},
            {b, innerHeight, .0f},
            {width, height, .0f},
            {width, innerHeight, .0f},
            {b, .0f, .0f},
            {width, .0f, .0f}
        }
    };

    Containers::Array<Vector3> normals(Containers::NoInit, 12);
    std::fill_n(normals.begin(), 12, Vector3::zAxis());

    Containers::Array<UnsignedInt> indices
        {
        Containers::InPlaceInit,
        {
                0u,1u,2u,
                0u,2u,3u,
                3u,5u,4u,
                3u,2u,5u,
                2u,6u,5u,
                2u,7u,6u,
                7u,8u,6u,
                7u,9u,8u,
                10u,9u,7u,
                10u,11u,9u
        }
        };

    Containers::Array<char> indexData{(char*)indices.release(), sizeof(UnsignedInt) * 10 * 3};

    Trade::MeshAttributeData vertexData{Trade::MeshAttribute::Position, VertexFormat::Vector3, vertices};
    Trade::MeshAttributeData normalData{Trade::MeshAttribute::Normal, VertexFormat::Vector3, normals};

    return MeshTools::interleave(
            Trade::MeshData(MeshPrimitive::Triangles, std::move(indexData), Trade::MeshIndexData{indices}, vertices.size()),
            {vertexData,normalData}
            );
}

void LoadPrimitives::load(ComboElement& element){
    switch (element.type) {
        case PrimitiveType::Capsule :
        {
            auto& optC = dynamic_cast<CapsuleOptions&>(*element.options);
            m_phasefieldData->original = Primitives::capsule3DSolid(optC.hemisphereRings, optC.cylinderRings, optC.segments, optC.segments);
            break;
        }
        case PrimitiveType::U :
        {
            auto& optU = dynamic_cast<UOptions&>(*element.options);
            m_phasefieldData->original = uShapedSquare(optU.width, optU.height, optU.innerWidth, optU.innerHeight);
            break;
        }
    }
    m_phasefieldData->status = PhasefieldData::Status::NewMesh;
}
