//
// Created by janos on 09.06.20.
//

#pragma once

#include "types.hpp"

class VisualizationProxy {

    void setFaceColors(ArrayView<double>& data);
    void setFaceColors(ArrayView<Mg::Color3ub>& data);
    void setVertexColors(ArrayView<double>& data);



    void* viewer;
};



