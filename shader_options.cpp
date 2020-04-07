//
// Created by janos on 07.04.20.
//

#include "shader_options.hpp"

#include <imgui.h>


void ShaderOptions::handleFlat() {

}


void ShaderOptions::handlePhong() {

}


void ShaderOptions::handleMeshVis() {

}

void ShaderOptions::drawImGui() {
    if (ImGui::TreeNode("Shader Options"))
    {
        const char* listbox_items[] = {
                "Phong Color Mapped",
                "Mesh Visualization",
                };

        static int listbox_item_current = 0;
        if(ImGui::ListBox("Shader Type\n", &listbox_item_current, listbox_items, IM_ARRAYSIZE(listbox_items), 4)){
            switch (listbox_item_current) {
                case 0 :
                    m_type = DrawableType::ColorMapPhong;
                    break;
                case 1 :
                    m_type = DrawableType::MeshVisualizer;
                    break;
                default:
                    break;
            }
            m_updateNode = true;
        }

        switch (listbox_item_current) {
            case 0 : handlePhong(); break;
            case 1 : handleFlat(); break;
            case 2 : handleMeshVis(); break;
            default : break;
        }
        ImGui::TreePop();
    }
}

void ShaderOptions::tickEvent(Scene &scene) {
    if(m_updateNode){
        m_phasefieldData.drawable = scene.setDrawableType("mesh", m_type);
        m_phasefieldData.type = m_type;
        m_updateNode = false;
    }
}
