//
// Created by janos on 09.06.20.
//

#pragma once

#include "Enums.h"

#include <Corrade/Containers/ArrayView.h>
#include <Corrade/Containers/Array.h>

namespace Phasefield {

struct VisualizationProxy {

    VisualizationProxy(class Viewer&);

    //void setFaceColors(Containers::ArrayView<double>& data);
    //void setFaceColors(Containers::ArrayView<Mg::Color3ub>& data);
    void setVertexColors(Containers::ArrayView<double>& data);

    void setVertexColors(class Tree& tree);

    void update();

    void setTag(Mg::Int tag);

    bool isActiveTag(Mg::Int tag) const;

    bool updateFaceTexture = false;
    bool updateVertexBuffer = false;
    //Containers::Array<Mg::Color3ub> faceColors;
    Viewer& viewer;

    Containers::Array<Mg::Color4> colors;

    Mg::Int activeTag = -1;
};

}