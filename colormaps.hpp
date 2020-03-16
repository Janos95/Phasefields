//
// Created by janos on 02.03.20.
//
#pragma once

#include <Magnum/Math/Color.h>
#include <Magnum/Magnum.h>
#include <Corrade/Containers/EnumSet.h>

#include <vector>
#include <algorithm>


class JetColorMap{
public:
     Magnum::Color4 operator()(double x);
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
