//
// Created by janos on 07.04.20.
//

#pragma once


#include "phasefield_data.hpp"
#include "viewer.hpp"

#include <Magnum/Math/Color.h>
#include <Magnum/Magnum.h>

enum class ColorMapType: Magnum::UnsignedShort {
    Turbo = 0,
    Magma = 1,
    Plasma = 2,
    Inferno = 3,
    Viridis = 4
};

using namespace Math::Literals;

std::unordered_map<ShaderType, Cr::Containers::Pointer<Mg::GL::AbstractShaderProgram>> makeShaders();
std::unordered_map<ColorMapType, Mg::GL::Texture2D> makeColorMapTextures();

void handleShaderOptions(
        Object3D* object,
        MeshDrawable*& drawable,
        DrawableGroup& drawableGroup,
        GL::Mesh& mesh,
        std::unordered_map<ShaderType, Containers::Pointer<GL::AbstractShaderProgram>>& shaders,
        std::unordered_map<ColorMapType, GL::Texture2D>& colorMapTextures,
        ColorMapType& colorMap);



