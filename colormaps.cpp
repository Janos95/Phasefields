//
// Created by janos on 02.03.20.
//

#include "colormaps.hpp"

#include <Magnum/Math/Color.h>
#include <Magnum/Magnum.h>

#include <algorithm> //thats ridiculous

using namespace Magnum;

namespace {

    float jet_colormap_red(float x) {
        if (x < 0.7f) {
            return 4.0f * x - 1.5f;
        } else {
            return -4.0f * x + 4.5f;
        }
    }

    float jet_colormap_green(float x) {
        if (x < 0.5f) {
            return 4.0f * x - 0.5f;
        } else {
            return -4.0f * x + 3.5f;
        }
    }

    float jet_colormap_blue(float x) {
        if (x < 0.3f) {
            return 4.0f * x + 0.5f;
        } else {
            return -4.0f * x + 2.5f;
        }
    }
}

Color4 JetColorMap::operator()(double x) {
    float r = std::clamp(jet_colormap_red(x), 0.0f, 1.0f);
    float g = std::clamp(jet_colormap_green(x), 0.0f, 1.0f);
    float b = std::clamp(jet_colormap_blue(x), 0.0f, 1.0f);
    return Color4(r, g, b, 1.0f);
}

