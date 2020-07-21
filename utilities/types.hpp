//
// Created by janos on 02.05.20.
//

#pragma once

#include <Magnum/Magnum.h>
#include <Magnum/SceneGraph/SceneGraph.h>
#include <Corrade/Containers/EnumSet.h>
#include <Corrade/Containers/Containers.h>

namespace Cr = Corrade;
namespace Mg = Magnum;
namespace Containers = Corrade::Containers;

using Object3D = Mg::SceneGraph::Object<Mg::SceneGraph::MatrixTransformation3D>;
using Scene3D = Mg::SceneGraph::Scene<Mg::SceneGraph::MatrixTransformation3D>;
using Drawable = Mg::SceneGraph::Drawable3D;
using DrawableGroup = Mg::SceneGraph::DrawableGroup3D;

enum class ColorMapType : Magnum::UnsignedInt {
    Turbo = 0,
    Magma = 1,
    Plasma = 2,
    Inferno = 3,
    Viridis = 4
};

enum class ShaderType : Magnum::UnsignedInt {
    FlatTextured = 0,
    PhongDiffuseColored = 1,
    PhongDiffuse = 2,
    MeshVisualizer = 3,
    MeshVisualizerPrimitiveId = 4,
    VertexColor = 8,
    Phong
};

enum class FunctionalType : Magnum::UnsignedInt {
    Undefined = 0,
    DirichletEnergy = 1,
    DoubleWellPotential = 2,
    Area = 3,
    Connectedness = 4,
};

enum class DrawableType : Magnum::Int {
    MeshVisualizer = 0,
    PhongDiffuse = 1,
    FlatTextured = 2,
    FaceColored = 3,
};



