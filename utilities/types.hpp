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

using Object3D = Mg::SceneGraph::Object<Mg::SceneGraph::MatrixTransformation3D>;
using Scene3D = Mg::SceneGraph::Scene<Mg::SceneGraph::MatrixTransformation3D>;
using Drawable = Mg::SceneGraph::Drawable3D;
using DrawableGroup = Mg::SceneGraph::DrawableGroup3D;

template<class T>
using View1D = Cr::Containers::StridedArrayView1D<T>;
template<class T>
using View2D = Cr::Containers::StridedArrayView2D<T>;
template<class T>
using View3D = Cr::Containers::StridedArrayView3D<T>;

using namespace Cr::Containers;

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
    Area1 = 3,
    Connectedness = 4,
    Area2 = 5,
};

enum class DrawableType : Magnum::Int {
    MeshVisualizer = 0,
    PhongDiffuse = 1,
    FlatTextured = 2,
    FaceColored = 3,
};

enum class VisualizationFlag : Magnum::UnsignedInt {
    Phasefield = 1u << 0u,
    GeodesicWeights = 1u << 1u,
    ConnectedComponents = 1u << 2u,
    Gradient = 1u << 3u,
};

using VisualizationFlags = Cr::Containers::EnumSet<VisualizationFlag>;
CORRADE_ENUMSET_OPERATORS(VisualizationFlags)


