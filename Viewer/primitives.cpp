//
// Created by janos on 04.04.20.
//

#include "primitives.hpp"
#include "Upload.h"
#include "custom_widgets.hpp"
//#include "polygonize_expression.hpp"

#include <Corrade/Containers/Pointer.h>
#include <Corrade/Containers/GrowableArray.h>

#include <Magnum/Primitives/Capsule.h>
#include <Magnum/MeshTools/Interleave.h>
#include <Magnum/MeshTools/Transform.h>
#include <Magnum/Math/Matrix4.h>

#include <imgui.h>
#include <misc/cpp/imgui_stdlib.h>


using namespace Corrade;
using namespace Magnum;

namespace {

struct ComboElement {
    std::string name;
    PrimitiveType type;
    Corrade::Containers::Pointer<AbstractPrimitiveOptions> options;
};

auto makeComboMap() {
    Containers::Array<ComboElement> map;
    Containers::arrayAppend(map, ComboElement{"U Shaped Square", PrimitiveType::U, Containers::pointer<UOptions>()});
    Containers::arrayAppend(map,
                            ComboElement{"Capsule", PrimitiveType::Capsule, Containers::pointer<CapsuleOptions>()});
    Containers::arrayAppend(map, ComboElement{"Implicit Function", PrimitiveType::ImplicitFunction,
                                              Containers::pointer<PolygonizationOptions>()});
    return map;
}

}


Trade::MeshData uShapedSquare(float width, float height, float innerWidth, float innerHeight) {
    float a = .5f*(width - innerWidth);
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
                            0u, 1u, 2u,
                            0u, 2u, 3u,
                            3u, 5u, 4u,
                            3u, 2u, 5u,
                            2u, 6u, 5u,
                            2u, 7u, 6u,
                            7u, 8u, 6u,
                            7u, 9u, 8u,
                            10u, 9u, 7u,
                            10u, 11u, 9u
                    }
            };

    auto size = indices.size();
    Trade::MeshIndexData indexData{indices};
    Containers::Array<char> data{reinterpret_cast<char*>(indices.release()), sizeof(UnsignedInt)*size};

    Trade::MeshAttributeData vertexData{Trade::MeshAttribute::Position, VertexFormat::Vector3, vertices};
    Trade::MeshAttributeData normalData{Trade::MeshAttribute::Normal, VertexFormat::Vector3, normals};

    return MeshTools::interleave(
            Trade::MeshData(MeshPrimitive::Triangles, std::move(data), indexData, vertices.size()),
            {vertexData, normalData}
    );
}

bool displayCapsuleOptions(CapsuleOptions& options) {
    constexpr static std::uint32_t step = 1;
    constexpr static float floatMax = 10.f;
    constexpr static float floatMin = .1f;

    constexpr static std::uint32_t ringsMin = 1u;
    constexpr static std::uint32_t segmentsMin = 3u;
    constexpr static std::uint32_t uintMax = 100u;
    bool hasChanged = false;
    hasChanged |= ImGui::SliderScalar("Numer of (face) rings for each heimsphere", ImGuiDataType_U32,
                                      &options.hemisphereRings, &ringsMin, &uintMax, "%u");
    hasChanged |= ImGui::SliderScalar("Number of (face) rings for cylinder", ImGuiDataType_U32, &options.cylinderRings,
                                      &ringsMin, &uintMax, "%u");
    hasChanged |= ImGui::SliderScalar("Number of (face) segments", ImGuiDataType_U32, &options.segments, &segmentsMin,
                                      &uintMax, "%u");
    hasChanged |= ImGui::SliderScalar("Radius of cylindinger", ImGuiDataType_Float, &options.radius, &floatMin,
                                      &floatMax, "%f");
    hasChanged |= ImGui::SliderScalar("Length of whole capsule", ImGuiDataType_Float, &options.length, &floatMin,
                                      &floatMax, "%f");

    return hasChanged;
}

bool displayUOptions(UOptions& options) {
    constexpr static std::uint32_t step = 1;
    constexpr static float floatMax = 10.f;
    constexpr static float floatMin = .1f;

    constexpr static std::uint32_t ringsMin = 1u;
    constexpr static std::uint32_t segmentsMin = 3u;
    constexpr static std::uint32_t uintMax = 100u;
    bool hasChanged = false;

    hasChanged |= ImGui::DragFloatRange2("Widths", &options.innerWidth, &options.width, .01f, .01f, 10.f, "Min: %.2f",
                                         "Max: %.2f");
    hasChanged |= ImGui::DragFloatRange2("Heights", &options.innerHeight, &options.height, .01f, .01f, 10.0f,
                                         "Min: %.2f", "Max: %.2f");
    return hasChanged;
}

void displayImplicitFunctionOptions(PolygonizationOptions& options, std::string& expression) {
    ImGui::InputTextMultiline("Implicit Function", &expression);
    constexpr static float minAngle = 10.f, maxAngle = 60.f;
    constexpr static float minBoundingRadius = 1., maxBoundingRadius = 10.f;
    constexpr static float minDistanceBound = 0.001f, maxDistanceBound = .5f;
    constexpr static float minRadiusBound = 0.001f, maxRadiusBound = .5f;
    ImGui::SliderScalar("Angle Bound", ImGuiDataType_Float, &options.angleBound, &minAngle, &maxAngle, "%f");
    ImGui::SliderScalar("Bounding Sphere Radius", ImGuiDataType_Float, &options.boundingSphereRadius,
                        &minBoundingRadius, &maxBoundingRadius, "%f");
    ImGui::SliderScalar("Distance Bound", ImGuiDataType_Float, &options.distanceBound, &minDistanceBound,
                        &maxDistanceBound, "%f");
    ImGui::SliderScalar("Radius Bound", ImGuiDataType_Float, &options.radiusBound, &minRadiusBound, &maxRadiusBound,
                        "%f");
}

bool handlePrimitive(Trade::MeshData& original, std::string& expression) {
    bool newMesh = false;

    static auto map = makeComboMap();
    static auto current = map.end();

    if(ImGui::BeginCombo("##combo", current != map.end() ? current->name.data() : nullptr)){
        for(auto it = map.begin(); it < map.end(); ++it){
            bool isSelected = (current ==
                               it); // You can store your selection however you want, outside or inside your objects
            if(ImGui::Selectable(it->name.c_str(), isSelected))
                current = it;
            if(isSelected)
                ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
        }
        ImGui::EndCombo();
    }

    auto buttonPressed = ImGui::Button("Load");
    ImGui::SameLine();
    static bool track = false;
    toggleButton("tracker", &track);

    bool hasChanged = false;
    if(current != map.end()){
        switch(current->type) {
            case PrimitiveType::Capsule : {
                auto& optC = dynamic_cast<CapsuleOptions&>(*current->options);
                hasChanged = displayCapsuleOptions(optC);
                if((hasChanged && track) || buttonPressed){
                    auto capsule = Primitives::capsule3DSolid(optC.hemisphereRings, optC.cylinderRings,
                                                              optC.segments, .5f*optC.length/optC.radius);
                    MeshTools::transformPointsInPlace(
                            Matrix4::scaling({optC.radius, optC.radius, optC.radius}),
                            capsule.mutableAttribute<Vector3>(Trade::MeshAttribute::Position));
                    original = std::move(capsule);
                    newMesh = true;
                }
                break;
            }
            case PrimitiveType::U : {
                auto& optU = dynamic_cast<UOptions&>(*current->options);
                hasChanged = displayUOptions(optU);
                if((hasChanged && track) || buttonPressed){
                    original = uShapedSquare(optU.width, optU.height, optU.innerWidth, optU.innerHeight);
                    newMesh = true;
                }
                break;
            }
            case PrimitiveType::ImplicitFunction :
                auto& optP = dynamic_cast<PolygonizationOptions&>(*current->options);
                displayImplicitFunctionOptions(optP, expression);
                if(buttonPressed){
                    //original = polygonizeExpression(expression, optP);
                    //newMesh = true;
                }
        }
        ImGui::Text("Track Options");
    }
    return newMesh;
}



//void LoadPrimitives::finalize(){
//    pd.meshData = preprocess(pd.original, CompileFlag::GenerateSmoothNormals|CompileFlag::AddTextureCoordinates);
//    auto vertexCount = pd.meshData.vertexCount();
//    Containers::arrayResize(pd.phasefield, vertexCount);
//    auto points = pd.meshData.attribute<Vector3>(Trade::MeshAttribute::Position);
//    Containers::arrayResize(pd.V, points.size());
//    for (int i = 0; i < vertexCount; ++i)
//        pd.V[i] = Vector3d(points[i]);
//    pd.F = pd.meshData.indicesAsArray();
//    Containers::arrayResize(pd.phasefield, vertexCount);
//    auto textureView = pd.meshData.mutableAttribute<Vector2>(Trade::MeshAttribute::TextureCoordinates);
//    for (int i = 0; i < vertexCount; ++i) {
//        textureView[i].x() = .5f * (pd.phasefield[i] + 1);
//    }
//    upload(pd, pd.meshData);
//    updateStatus = true;
//}
