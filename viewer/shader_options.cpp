//
// Created by janos on 07.04.20.
//


#include "shader_options.hpp"
#include "custom_widgets.hpp"

#include <Corrade/Containers/Array.h>
#include <Corrade/Containers/Pointer.h>

#include <Magnum/DebugTools/ColorMap.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/ImageView.h>
#include <Magnum/GL/Renderer.h>

#include <imgui.h>

using namespace Corrade;
using namespace Magnum;

namespace {
    struct ComboElement {
        std::string name;
        ColorMapType type;
    };
    auto makeComboElementMap(){
        Containers::Array<ComboElement> map{Containers::InPlaceInit, {
                {"Turbo", ColorMapType::Turbo},
                {"Magma", ColorMapType::Magma},
                {"Plasma", ColorMapType::Plasma},
                {"Inferno", ColorMapType::Inferno},
                {"Viridis", ColorMapType::Viridis}
        }};
        return map;
    }
}


void drawMeshVisualizerOptions(MeshVisualizerDrawable& drawable) {
    constexpr Float min = 0.f, max = 10.f;
    ImGui::DragScalar("Wireframe Width", ImGuiDataType_Float, &drawable.wireframeWidth, .01f, &min, &max, "%f", 1);
    ImGui::ColorEdit3("Wireframe Color", drawable.wireframeColor.data());
    ImGui::ColorEdit3("Color", drawable.color.data());
    ImGui::DragScalar("Smoothness", ImGuiDataType_Float, &drawable.smoothness, .01f, &min, &max, "%f", 1);
}

bool drawColorMapOptions(ColorMapType& mapType) {
    static auto map = makeComboElementMap();
    auto current = std::find_if(map.begin(), map.end(), [&](auto& el){ return el.type == mapType; });
    bool newMap = false;
    if (ImGui::BeginCombo("##combo", current->name.data())) {
        for (auto it = map.begin(); it < map.end(); ++it) {
            bool isSelected = (current == it);
            if (ImGui::Selectable(it->name.c_str(), isSelected)) {
                current = it;
                mapType = it->type;
                newMap = true;
            }

            if (isSelected)
                ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
    return newMap;
}


std::unordered_map<ShaderType, Cr::Containers::Pointer<Mg::GL::AbstractShaderProgram>> makeShaders() {
    std::unordered_map<ShaderType, Cr::Containers::Pointer<Mg::GL::AbstractShaderProgram>> map;

    auto Wireframe = Shaders::MeshVisualizer3D::Flag::Wireframe;
    auto PrimitiveColored = Shaders::MeshVisualizer3D::Flag::PrimitiveId;

    map.emplace(ShaderType::PhongDiffuse, new Shaders::Phong(Shaders::Phong::Flag::DiffuseTexture));
    map.emplace(ShaderType::MeshVisualizer, new Shaders::MeshVisualizer3D(Wireframe));
    map.emplace(ShaderType::MeshVisualizerObjectId, new Shaders::MeshVisualizer3D(PrimitiveColored));
    map.emplace(ShaderType::FlatTextured, new Shaders::Flat3D(Shaders::Flat3D::Flag::Textured));
    map.emplace(ShaderType::VertexColor, new Shaders::VertexColor3D());

    return map;
}


std::unordered_map<ColorMapType, Mg::GL::Texture2D> makeColorMapTextures(){
    std::unordered_map<ColorMapType, Mg::GL::Texture2D> map;
    using L = std::initializer_list<std::pair<ColorMapType, Containers::StaticArrayView<256, const Vector3ub>>>;
    for(auto&& [type, colorMap] : L{
            {ColorMapType::Turbo, Magnum::DebugTools::ColorMap::turbo()},
            {ColorMapType::Magma, Magnum::DebugTools::ColorMap::magma()},
            {ColorMapType::Plasma, Magnum::DebugTools::ColorMap::plasma()},
            {ColorMapType::Inferno, Magnum::DebugTools::ColorMap::inferno()},
            {ColorMapType::Viridis, Magnum::DebugTools::ColorMap::viridis()}
    })
    {
        const Magnum::Vector2i size{Magnum::Int(colorMap.size()), 1};
        GL::Texture2D texture;
        texture.setMinificationFilter(Magnum::SamplerFilter::Linear)
                .setMagnificationFilter(Magnum::SamplerFilter::Linear)
                .setWrapping(Magnum::SamplerWrapping::ClampToEdge) // or Repeat
                .setStorage(1, Magnum::GL::TextureFormat::RGB8, size) // or SRGB8
                .setSubImage(0, {}, ImageView2D{Magnum::PixelFormat::RGB8Srgb, size, colorMap});
        map.emplace(type, std::move(texture));
    }
    return map;
}


