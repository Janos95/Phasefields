//
// Created by janos on 27.02.20.
//

#pragma once

#include "object.hpp"

#include <Magnum/Math/Matrix4.h>
#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Shaders/Flat.h>

class DefaultCallback
{
public:

    explicit DefaultCallback(Object&, Magnum::Shaders::Flat3D&);

    void operator()(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera);

    DefaultCallback(const DefaultCallback&) = default;

private:

    Magnum::GL::Mesh& m_mesh;
    Magnum::GL::Texture2D* m_texture = nullptr;
    Magnum::Color4 m_color;
    Magnum::Shaders::Flat3D& m_shader;
};



