//
// Created by janos on 07.04.20.
//

#include "shader_options.hpp"

#include <imgui.h>


bool ShaderOptions::displayPhongOptions() {

}


bool ShaderOptions::displayFlatOptions() {

}


bool ShaderOptions::displayMeshVisualizerOptions() {

}

ShaderOptions::ShaderOptions(PhasefieldData &data):
    m_phasefieldData(data)
{
    m_shaders.emplace_back(new Shaders::Flat3D{});
    m_shaders.emplace_back(new Shaders::Flat3D{Shaders::Flat3D::Flag::Textured});
    m_shaders.emplace_back(new Shaders::VertexColor3D{});
    m_shaders.emplace_back(new Shaders::MeshVisualizer3D{Shaders::MeshVisualizer3D::Flag::Wireframe});
    m_shaders.emplace_back(new Shaders::Phong{Shaders::Phong::Flag::DiffuseTexture, 2});
}

void ShaderOptions::drawImGui(Viewer& viewer) {
    if (ImGui::TreeNode("Shader Options"))
    {
        const char* listbox_items[] = {
                "Mesh Visualization",
                "Phong Color Mapped",
                "Flat Color Mapped"
                };

        static DrawableType type = DrawableType::ColorMapPhong;
        bool newShader = ImGui::ListBox("Shader Type\n", (int*)&type, listbox_items, IM_ARRAYSIZE(listbox_items), 4);
        if(newShader)
            m_node->features().clear();

        switch (type) {
            case DrawableType::MeshVisualizer :
                if(displayPhongOptions() || newShader)
                    new MeshVisualizerDrawable(*m_node, *m_shaders[4], &viewer.drawableGroup);
            case DrawableType::ColorMapPhong :
                if(displayFlatOptions() || newShader)
                    new ColorMapPhongDrawable(*m_node, *m_shaders[4], &viewer.drawableGroup);
            case DrawableType::ColorMapFlat :
                if(displayMeshVisualizerOptions() || newShader)
                    new ColorMapFlatDrawable(*m_node, *m_shaders[4], &viewer.drawableGroup);
            default : break;
        }

        ImGui::TreePop();
    }
}

