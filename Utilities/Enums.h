//
// Created by janos on 02.05.20.
//

#pragma once

#include "SmartEnum.h"

#include <Corrade/Containers/Containers.h>
#include <Magnum/Types.h>

namespace Cr = Corrade;
namespace Mg = Magnum;

namespace Containers = Corrade::Containers;

SMART_ENUM(ColorMapType, Mg::UnsignedInt,
        Turbo,
        Magma,
        Plasma,
        Inferno,
        Viridis)

enum class ShaderType : Magnum::UnsignedInt {
    FlatTextured = 0,
    PhongDiffuseColored = 1,
    PhongDiffuse = 2,
    MeshVisualizer = 3,
    MeshVisualizerPrimitiveId = 4,
    VertexColor = 8,
    Phong
};

SMART_ENUM(FunctionalType, Mg::UnsignedInt,
           DirichletEnergy,
           DoubleWellPotential,
           AreaRegularizer,
           ConnectednessConstraint,
           DiffuseYamabe,
           HierarchicalRegularization,
           Unknown)

enum class DrawableType : Magnum::Int {
    MeshVisualizer = 0,
    PhongDiffuse = 1,
    FlatTextured = 2,
    FaceColored = 3,
};



