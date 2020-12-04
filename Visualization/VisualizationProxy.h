//
// Created by janos on 09.06.20.
//

#pragma once

#include "Enums.h"
#include "UniqueFunction.h"
#include "Surface.h"


//TODO UniqueFunction blacks out with reference to incomplete type ??
#include "Tree.h"

namespace Phasefield {

SMART_ENUM(VisOption, size_t,
           Phasefield,
           Weight,
           Segmentation)

/**
 * @param n numbers of colors you want
 * @return array of hopefully visually distinct colors
 */
Array<Color4>& getColors(size_t n = Invalid);

void optimizeColors(Array<size_t> const& neighbors, Array<size_t> const& starts);

/**
 * The main purpose of this class is to avoid including the whole Viewer.h
 * everywhere we want to manipulate some of the drawing state.
 */
struct VisualizationProxy {

    enum class ShaderConfig {
        VertexColors,
        ColorMaps
    };

    ShaderConfig shaderConfig = ShaderConfig::ColorMaps;

    VisualizationProxy(class Viewer&);

    //void setFaceColors(Containers::ArrayView<double>& data);
    //void setFaceColors(Containers::ArrayView<Mg::Color3ub>& data);

    void redraw();

    void upload();

    bool updateVertexBuffer = false;

    //Containers::Array<Mg::Color3ub> faceColors;
    Viewer& viewer;

    void drawSegmentation();

    void drawValues(VertexDataView<double> const& values);

    void drawValuesNormalized(VertexDataView<double> const& values);

    void setDefaultCallback();

    void setCallbacks(UniqueFunction<void(Node)>, UniqueFunction<void()>);

    UniqueFunction<void(Node)> drawCb;
    UniqueFunction<void()> releaseCb;

    VisOption::Value option = VisOption::Phasefield;

    bool customMapping = false;
    double scale = 1;
    double offset = 0;
    int level = 0;
};

}