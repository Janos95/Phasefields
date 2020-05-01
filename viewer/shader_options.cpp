//
// Created by janos on 07.04.20.
//

#include "shader_options.hpp"
#include "custom_widgets.hpp"

#include <Magnum/DebugTools/ColorMap.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/ImageView.h>
#include <Magnum/GL/Renderer.h>

#include <imgui.h>


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


bool drawMeshVisualizerOptions(MeshVisualizerDrawable& drawable) {

    bool updateShader = false;
    ImGui::Text("Draw Normals");
    ImGui::SameLine();
    updateShader |= toggleButton("normals", &drawable.drawNormals);
    ImGui::SameLine();
    ImGui::Text("Draw Tangent Space");
    ImGui::SameLine();
    updateShader |= toggleButton("tangentspace", &drawable.drawTangentSpace);
    if(drawable.drawNormals || drawable.drawTangentSpace)
        Debug{} << "drawing normals or tangent space";

    auto& shader = *drawable.shader;
    constexpr Float min = 0.f, max = 10.f;

    ImGui::DragScalar("Wireframe Width", ImGuiDataType_Float, &drawable.wireframeWidth, .01f, &min, &max, "%f", 1);
    ImGui::DragScalar("Line Width", ImGuiDataType_Float, &drawable.lineWidth, .01f, &min, &max, "%f", 1);
    ImGui::DragScalar("Line Length", ImGuiDataType_Float, &drawable.lineLength, .01f, &min, &max, "%f", 1);
    ImGui::ColorEdit3("Wireframe Color", drawable.wireframeColor.data());
    ImGui::ColorEdit3("Color", drawable.color.data());
    ImGui::DragScalar("Smoothness", ImGuiDataType_Float, &drawable.smoothness, .01f, &min, &max, "%f", 1);

    return updateShader;
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
    auto Normal = Shaders::MeshVisualizer3D::Flag::NormalDirection | Wireframe;
    auto Bitangents = Shaders::MeshVisualizer3D::Flag::BitangentFromTangentDirection | Wireframe;
    auto Tangents = Shaders::MeshVisualizer3D::Flag::TangentDirection | Wireframe;
    auto ObjectId = Shaders::MeshVisualizer3D::Flag::InstancedObjectId | Wireframe;

    map.emplace(ShaderType::PhongDiffuse, new Shaders::Phong(Shaders::Phong::Flag::DiffuseTexture));
    map.emplace(ShaderType::MeshVisualizer, new Shaders::MeshVisualizer3D(Wireframe));
    map.emplace(ShaderType::MeshVisualizerNormal, new Shaders::MeshVisualizer3D(Normal));
    map.emplace(ShaderType::MeshVisualizerTangent, new Shaders::MeshVisualizer3D(Tangents | Bitangents));
    map.emplace(ShaderType::MeshVisualizerFull, new Shaders::MeshVisualizer3D(Tangents | Bitangents | Normal));
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

void handleShaderOptions(
        Object3D* object,
        MeshDrawable*& drawable,
        DrawableGroup& drawableGroup,
        GL::Mesh& mesh,
        std::unordered_map<ShaderType, Containers::Pointer<GL::AbstractShaderProgram>>& shaders,
        std::unordered_map<ColorMapType, GL::Texture2D>& colorMapTextures,
        ColorMapType& colorMap){

    static constexpr std::pair<char const*, DrawableType> drawablesMap[] = {
            {"Mesh Visualization", DrawableType::MeshVisualizer},
            {"Phong Diffuse Textured", DrawableType::PhongDiffuse},
            {"Flat", DrawableType::FlatTextured}
    };

    static auto* current = drawablesMap + 1;
    bool newDrawable = false;
    if (ImGui::BeginCombo("##combo", current->first)) {
        for (auto it = drawablesMap; it != drawablesMap + 3; ++it) {
            bool isSelected = (current == it);
            if (ImGui::Selectable(it->first, isSelected)){
                newDrawable = true;
                current = it;
            }
            if (isSelected)
                ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
        }
        ImGui::EndCombo();
    }

    if(newDrawable)
        object->features().clear();

    switch (current->second) {
        case DrawableType::MeshVisualizer : {
            if (newDrawable)
                drawable = new MeshVisualizerDrawable(*object, mesh, &drawableGroup);
            auto &meshVis = dynamic_cast<MeshVisualizerDrawable &>(*drawable);
            if(drawMeshVisualizerOptions(meshVis) || newDrawable){
                if(meshVis.drawTangentSpace && meshVis.drawNormals){
                    meshVis.shader = dynamic_cast<Shaders::MeshVisualizer3D*>(shaders[ShaderType::MeshVisualizerFull].get());
                } else if(meshVis.drawTangentSpace){
                    meshVis.shader = dynamic_cast<Shaders::MeshVisualizer3D*>(shaders[ShaderType::MeshVisualizerTangent].get());
                } else if(meshVis.drawNormals){
                    meshVis.shader = dynamic_cast<Shaders::MeshVisualizer3D*>(shaders[ShaderType::MeshVisualizerNormal].get());
                } else {
                    meshVis.shader = dynamic_cast<Shaders::MeshVisualizer3D*>(shaders[ShaderType::MeshVisualizer].get());
                }
            }
            break;
        }
        case DrawableType::PhongDiffuse : {
            if (newDrawable)
                drawable = new PhongDiffuseDrawable(*object, mesh, *shaders[ShaderType::PhongDiffuse], &drawableGroup);
            auto& phongVis = dynamic_cast<PhongDiffuseDrawable&>(*drawable);
            if(drawColorMapOptions(colorMap) || newDrawable)
                phongVis.texture = &colorMapTextures[colorMap];
            break;
        }
        case DrawableType::FlatTextured : {
            if (newDrawable)
                drawable = new FlatTexturedDrawable(*object, mesh, *shaders[ShaderType::FlatTextured], &drawableGroup);
            auto& flatVis = dynamic_cast<FlatTexturedDrawable&>(*drawable);
            if(drawColorMapOptions(colorMap) || newDrawable)
                flatVis.texture = &colorMapTextures[colorMap];
            break;
        }
    }
}

