//
// Created by janos on 08.05.20.
//

#pragma once

#include "Enums.h"

#include <Magnum/Primitives/Cylinder.h>
#include <Magnum/Math/Matrix4.h>
#include <Magnum/Math/Matrix3.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/GL/Buffer.h>

#include <Corrade/Containers/Array.h>


namespace Mg = Magnum;
namespace Cr = Corrade;

using namespace Mg::Math::Literals;

struct InstanceData {
    Mg::Matrix4 tf;
    Mg::Matrix3 normalMatrix;
    Mg::Color3 color;
};

/**
 * @todo Instead of instancing a better options would be to simply upload
 * @todo all edges as a line set and write a custom shader which for each
 * @todo edge generates a sreen aligned quad in a geometry shader and some normals
 * @todo so that we can than add some plasticity. For junctions either use blending
 * @todo or compute an sdf and use a second render pass to decide which normals
 * @todo to use for shading.
 */

struct Paths {
    explicit Paths();
    void draw(Mg::Matrix4 const& transformation, Mg::Matrix4 const& projection);

    Mg::Shaders::Phong shader;
    Mg::GL::Mesh cylinder;
    Mg::GL::Buffer instanceBuffer;
    Mg::Containers::Array<InstanceData> instanceData;
    Mg::Containers::Array<InstanceData> instanceDataTransformed;
    bool drawPaths = false;
};


