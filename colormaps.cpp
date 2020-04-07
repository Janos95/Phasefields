//
// Created by janos on 02.03.20.
//

#include "colormaps.hpp"

#include <Magnum/Math/Color.h>
#include <Magnum/Magnum.h>
#include <Magnum/DebugTools/ColorMap.h>
#include <Magnum/GL/TextureFormat.h>
#include <Magnum/PixelFormat.h>
#include <Magnum/ImageView.h>


#include <imgui.h>

using namespace Magnum;


ColorMap::ColorMap(PhasefieldData &data):
    phasefieldData(data),
    map{
            {"Turbo", ColorMapType::Turbo},
            {"Magma", ColorMapType::Magma},
            {"Plasma", ColorMapType::Plasma},
            {"Inferno", ColorMapType::Inferno},
            {"Viridis", ColorMapType::Viridis}
    }
{
    current = map.begin();
}

void ColorMap::drawImGui(){
    if (ImGui::TreeNode("Color Maps"))
    {
        if (ImGui::BeginCombo("##combo", current->name.data())) {
            for (auto it = map.begin(); it < map.end(); ++it){
                bool isSelected = (current == it); // You can store your selection however you want, outside or inside your objects
                if (ImGui::Selectable(it->name.c_str(), isSelected)){
                    auto drawable = phasefieldData.drawable;
                    auto type = phasefieldData.type;
                    if(drawable && it != current && type == DrawableType::ColorMapPhong){
                        current = it;
                        switch(phasefieldData.type){
                            case DrawableType::ColorMapPhong :
                                dynamic_cast<ColorMapPhongDrawable*>(drawable)->setColorMapping(current->type);
                                break;
                            default:
                                break;
                        }
                    }
                }

                if (isSelected)
                    ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
            }
            //@todo imgui display image of colormap
            ImGui::EndCombo();
        }

        ImGui::TreePop();
    }
}


