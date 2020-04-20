//
// Created by janos on 02.03.20.
//
#pragma once

#include "phasefield_data.hpp"
#include "viewer.hpp"

#include <Corrade/Containers/EnumSet.h>

#include <Magnum/Math/Color.h>
#include <Magnum/Magnum.h>
#include <Magnum/DebugTools/ColorMap.h>

#include <vector>
#include <algorithm>




struct ColorMap : Viewer::AbstractEventHandler {

    struct ComboElement{
        ComboElement(std::string n, ColorMapType t) : name(std::move(n)), type(t) {}
        std::string name;
        ColorMapType type;
    };

    explicit ColorMap(PhasefieldData& data);

    PhasefieldData& phasefieldData;

    std::vector<ComboElement> map;
    std::vector<ComboElement>::iterator current;
    bool updated = true;

    void drawImGui(Viewer&) override;
};

class SegmentationColorMap{
public:
    Magnum::Color4 operator()(double x){
        auto it = std::find_if(interfaces.begin(), interfaces.end(), InInterface{x});
        if(it != interfaces.end())
            return it->c;
        return undefined;
    }
    
    struct Interface{
        double a, b;
        Magnum::Color4 c;
    };
    
    struct InInterface{
        double x;
        bool operator()(Interface const& interf){
            return interf.a <= x && x <= interf.b;
        }
    };
   
    static constexpr Magnum::Color4 undefined = Magnum::Color4(1,1,1,1);
    std::vector<Interface> interfaces;
};
