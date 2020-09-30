//
// Created by janos on 09.06.20.
//

#pragma once

#include "Enums.h"
#include "Mesh.h"
#include "UniqueFunction.h"

#include <Corrade/Containers/Array.h>

namespace Phasefield {

SMART_ENUM(VisOption, size_t,
           Phasefield,
           Weight,
           Segmentation)

/**
 * @param n numbers of colors you want
 * @return array of hopefully visually distinct colors
 */
Array<Color4>& getColors(size_t n);

struct VisualizationProxy {

    VisualizationProxy(class Viewer&);

    //void setFaceColors(Containers::ArrayView<double>& data);
    //void setFaceColors(Containers::ArrayView<Mg::Color3ub>& data);

    void redraw();

    void upload();

    bool updateVertexBuffer = false;

    //Containers::Array<Mg::Color3ub> faceColors;
    Viewer& viewer;

    void drawSegmentation();

    void drawValues(VertexDataView<double> values, UniqueFunction<double(double)> tf);

    void setDefaultCallback();

    void setCallbacks(UniqueFunction<void(Viewer*)>, UniqueFunction<void()>);

    UniqueFunction<void(Viewer*)> drawCb;
    UniqueFunction<void()> releaseCb;

    VisOption::Value option;
    bool isDefaultCallback = true;

    using P = std::pair<const char*, Array<Color4>>;
    Array<P> colorMaps;
};

}