//
// Created by janos on 08.05.20.
//

#pragma once

#include "types.hpp"
#include "drawables.hpp"

#include <Magnum/Primitives/Cylinder.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/Matrix3.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Trade/MeshData.h>
#include <Corrade/Containers/GrowableArray.h>

using namespace Mg::Math::Literals;

namespace Mg = Magnum;
namespace Cr = Corrade;

struct InstanceData {
    Mg::Matrix4 tf;
    Mg::Matrix3 normalMatrix;
    Mg::Color3 color;
};

struct Paths : Object3D, Drawable {
    explicit Paths(Object3D* parent, Mg::SceneGraph::DrawableGroup3D& drawables);
    void draw(const Mg::Matrix4& transformation, Mg::SceneGraph::Camera3D& camera) override;

    Mg::Shaders::Phong shader;
    Mg::GL::Mesh cylinder;
    Mg::GL::Buffer instanceBuffer;
    Mg::Containers::Array<InstanceData> instanceData;
    Mg::Containers::Array<InstanceData> instanceDataTransformed;
    bool drawPaths = false;
};


