//
// Created by janos on 8/26/20.
//

#pragma once

#include <Magnum/Magnum.h>
#include <Magnum/Shaders/Shaders.h>
#include <Magnum/GL/GL.h>

using namespace Corrade::Containers;
using namespace Magnum::Shaders;

using Magnum::UnsignedInt;
using Magnum::Float;
using Magnum::Double;
using Magnum::Int;

using Magnum::Vector3;
using Magnum::Vector2;
using Magnum::Vector2i;
using Magnum::Vector2ui;
using Magnum::Vector3ub;
using Magnum::Vector3ui;
using Magnum::Vector3d;
using Magnum::Matrix3;
using Magnum::Matrix3d;
using Magnum::Matrix4;
using Magnum::Color3;
using Magnum::Color3ub;
using Magnum::Color4;
using Magnum::Color4ub;
using Magnum::Range2D;
using Magnum::Range2Di;

using Magnum::Rad;
using Magnum::Radd;
using Magnum::Deg;

using Magnum::Debug;

namespace Math = Magnum::Math;
namespace GL = Magnum::GL;

namespace Phasefield {
constexpr size_t Invalid = ~size_t{0};
}
