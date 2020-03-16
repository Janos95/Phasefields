//
// Created by janos on 02.02.20.
//

#pragma once

#include <Corrade/Containers/Optional.h>

#include <Magnum/GL/Texture.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>

#include <memory>

struct Object
{
    Magnum::GL::Buffer vertices, indices;
    Magnum::GL::Mesh mesh;
    std::unique_ptr<Magnum::GL::Texture2D> textureDiffuse = nullptr;
    std::unique_ptr<Magnum::GL::Texture2D> textureSpecular = nullptr;
    Magnum::Color4 color = Magnum::Color4::blue();
};
