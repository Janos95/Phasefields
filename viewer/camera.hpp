//
// Created by janos on 02.02.20.
//

#pragma once

#include <Magnum/Math/Matrix4.h>

struct Camera
{
    Magnum::Matrix4 transformation;
    Magnum::Matrix4 projection;
};

