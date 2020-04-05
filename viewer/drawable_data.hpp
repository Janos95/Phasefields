//
// Created by janos on 02.02.20.
//

#pragma once

#include <Corrade/Containers/Optional.h>

#include <Magnum/GL/Texture.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Trade/MeshData.h>

#include <memory>

struct DrawableData
{
    Magnum::GL::Buffer vertices, indices;
    Magnum::GL::Mesh mesh;
    Magnum::Trade::MeshData meshData{Magnum::MeshPrimitive::Points, 0};
};

struct UniformColorDrawableData : virtual DrawableData
{
    Magnum::Color4 color;
};

struct TexturedDrawableData : virtual DrawableData
{
    Magnum::GL::Texture2D texture;
};






