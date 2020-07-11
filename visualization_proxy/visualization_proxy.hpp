//
// Created by janos on 09.06.20.
//

#pragma once

#include "../utilities/types.hpp"

#include <Corrade/Containers/ArrayView.h>
#include <Corrade/Containers/Array.h>

class Viewer;

struct VisualizationProxy {

    VisualizationProxy(Viewer&);

    void setFaceColors(Containers::ArrayView<double>& data);
    void setFaceColors(Containers::ArrayView<Mg::Color3ub>& data);
    void setVertexColors(Containers::ArrayView<double>& data);
    void resize();

    void update();

    bool updateFaceTexture = false;
    bool updateVertexUvs = false;
    Containers::Array<Mg::Color3ub> faceColors;
    Viewer& viewer;
};



