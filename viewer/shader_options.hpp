//
// Created by janos on 07.04.20.
//

#pragma once

#include "types.hpp"
#include "drawables.hpp"

#include <Magnum/Math/Color.h>
#include <Magnum/Magnum.h>
#include <Magnum/GL/Mesh.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

using namespace Mg::Math::Literals;

std::unordered_map<ShaderType, Cr::Containers::Pointer<Mg::GL::AbstractShaderProgram>> makeShaders();
std::unordered_map<ColorMapType, Mg::GL::Texture2D> makeColorMapTextures();

void drawMeshVisualizerOptions(MeshVisualizerDrawable& drawable);

bool drawColorMapOptions(ColorMapType& mapType);
