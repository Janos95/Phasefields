//
// Created by janos on 02.02.20.
//

#pragma once

#include "scene_graph_node.hpp"

#include <Corrade/Containers/Optional.h>

#include <Magnum/GL/Texture.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Magnum.h>
#include <Magnum/Math/Color.h>

struct Object
{
    Magnum::GL::Buffer vertices, indices;
    Magnum::GL::Mesh mesh;
    Corrade::Containers::Optional<Magnum::GL::Texture2D> texture;
    Magnum::Color4 color = Magnum::Color4{1};
    SceneGraphNode* node;
};
