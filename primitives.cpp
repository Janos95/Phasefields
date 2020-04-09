//
// Created by janos on 04.04.20.
//

#include "primitives.hpp"
#include "upload.hpp"
#include "toggle_button.hpp"

#include <Corrade/Utility/Algorithms.h>

#include <Magnum/Primitives/Capsule.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/Transform.h>
#include <Magnum/MeshTools/RemoveDuplicates.h>
#include <Magnum/Math/Matrix4.h>

#include <imgui.h>

using ComboElement = LoadPrimitives::ComboElement;

namespace {
    struct CompareByName{
        bool operator()(ComboElement const& e1, ComboElement const& e2){
            return e1.name < e2.name;
        }
    };

    auto makeComboMap(){
        std::vector<ComboElement> map;
        map.push_back(ComboElement{"U Shaped Square", PrimitiveType::U, std::make_unique<LoadPrimitives::UOptions>()});
        map.push_back(ComboElement{"Capsule", PrimitiveType::Capsule, std::make_unique<LoadPrimitives::CapsuleOptions>()});
        std::sort(map.begin(), map.end(), CompareByName{});
        return map;
    }

}

bool displayCapsuleOptions(LoadPrimitives::CapsuleOptions& options){
    constexpr static std::uint32_t step = 1;
    constexpr static float floatMax = 10.f;
    constexpr static float floatMin = .1f;

    constexpr static std::uint32_t ringsMin = 1u;
    constexpr static std::uint32_t segmentsMin = 3u;
    constexpr static std::uint32_t uintMax = 100u;
    bool hasChanged = false;
    hasChanged |= ImGui::SliderScalar("Numer of (face) rings for each heimsphere", ImGuiDataType_U32, &options.hemisphereRings, &ringsMin, &uintMax, "%u");
    hasChanged |= ImGui::SliderScalar("Number of (face) rings for cylinder", ImGuiDataType_U32, &options.cylinderRings, &ringsMin, &uintMax, "%u");
    hasChanged |= ImGui::SliderScalar("Number of (face) segments", ImGuiDataType_U32, &options.segments, &segmentsMin, &uintMax, "%u");
    hasChanged |= ImGui::SliderScalar("Radius of cylindinger", ImGuiDataType_Float, &options.radius, &floatMin, &floatMax, "%f");
    hasChanged |= ImGui::SliderScalar("Length of whole capsule", ImGuiDataType_Float, &options.length, &floatMin, &floatMax, "%f");
    return hasChanged;
}

bool displayUOptions(LoadPrimitives::UOptions& options){
    constexpr static std::uint32_t step = 1;
    constexpr static float floatMax = 10.f;
    constexpr static float floatMin = .1f;

    constexpr static std::uint32_t ringsMin = 1u;
    constexpr static std::uint32_t segmentsMin = 3u;
    constexpr static std::uint32_t uintMax = 100u;
    bool hasChanged = false;

    hasChanged |= ImGui::DragFloatRange2("Widths", &options.innerWidth, &options.width, .01f, .01f, 10.f, "Min: %.2f", "Max: %.2f");
    hasChanged |= ImGui::DragFloatRange2("Heights", &options.innerHeight, &options.height, .01f, .01f, 10.0f, "Min: %.2f", "Max: %.2f");
    return hasChanged;
}

void LoadPrimitives::drawImGui(){
    if (ImGui::TreeNode("Primitives"))
    {
        static auto map = makeComboMap();
        static auto current = map.end();

        if (ImGui::BeginCombo("##combo", current != map.end() ? current->name.data() : nullptr)) {
            for (auto it = map.begin(); it < map.end(); ++it){
                bool isSelected = (current == it); // You can store your selection however you want, outside or inside your objects
                if (ImGui::Selectable(it->name.c_str(), isSelected))
                    current = it;
                if (isSelected)
                    ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
            }

            ImGui::EndCombo();
        }

        bool hasChanged = false;
        if(current != map.end()){
            switch (current->type) {
                case PrimitiveType::Capsule :
                    hasChanged = displayCapsuleOptions(dynamic_cast<CapsuleOptions&>(*current->options));
                    break;
                case PrimitiveType::U :
                    hasChanged = displayUOptions(dynamic_cast<UOptions&>(*current->options));
                    break;
                default:
                    break;
            }

            toggleButton("Track Options,", &track);
        }

        if(((hasChanged && track) || ImGui::Button("Load")) && current != map.end()){
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

    auto size = indices.size();
    Trade::MeshIndexData indexData{indices};
    Containers::Array<char> data{reinterpret_cast<char*>(indices.release()), sizeof(UnsignedInt) * size};

    Trade::MeshAttributeData vertexData{Trade::MeshAttribute::Position, VertexFormat::Vector3, vertices};
    Trade::MeshAttributeData normalData{Trade::MeshAttribute::Normal, VertexFormat::Vector3, normals};

    return MeshTools::interleave(
            Trade::MeshData(MeshPrimitive::Triangles, std::move(data), indexData, vertices.size()),
            {vertexData,normalData}
            );
}

void LoadPrimitives::load(ComboElement& element){
    switch (element.type) {
        case PrimitiveType::Capsule :
        {
            auto& optC = dynamic_cast<CapsuleOptions&>(*element.options);
            auto capsule = Primitives::capsule3DSolid(optC.hemisphereRings, optC.cylinderRings, optC.segments, .5f * optC.length/optC.radius, Primitives::CapsuleFlag::Tangents);
            MeshTools::transformPointsInPlace(Math::Matrix4<float>::scaling({optC.radius,optC.radius,optC.radius}), capsule.mutableAttribute<Vector3>(Trade::MeshAttribute::Position));
            phasefieldData.original = std::move(capsule);
            break;
        }
        case PrimitiveType::U :
        {
            auto& optU = dynamic_cast<UOptions&>(*element.options);
            phasefieldData.original = uShapedSquare(optU.width, optU.height, optU.innerWidth, optU.innerHeight);
            break;
        }
    }
    phasefieldData.meshData = preprocess(phasefieldData.original, CompileFlag::GenerateSmoothNormals|CompileFlag::AddTextureCoordinates);
    phasefieldData.V = phasefieldData.meshData.positions3DAsArray();
    phasefieldData.F = phasefieldData.meshData.indicesAsArray();
    Containers::arrayResize(phasefieldData.phasefield, phasefieldData.V.size());
    phasefieldData.status = PhasefieldData::Status::NewMesh;
}
