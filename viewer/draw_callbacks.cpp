//
// Created by janos on 27.02.20.
//

#include "draw_callbacks.hpp"
#include "object.hpp"

#include <Magnum/Shaders/Flat.h>

using namespace Magnum;
using namespace Corrade;

FlatCallback::FlatCallback(Object& object, Shaders::Flat3D& shader) :
        m_mesh(object.mesh),
        m_color(object.color),
        m_shader(shader)
{
    if(m_shader.flags() & Shaders::Flat3D::Flag::Textured){
        CORRADE_INTERNAL_ASSERT(object.textureDiffuse);
        m_texture = object.textureDiffuse.get();
    }
}

void FlatCallback::operator()(const Matrix4& transformationMatrix, SceneGraph::Camera3D& camera){
    m_shader.setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix);

    if (m_shader.flags() & Shaders::Flat3D::Flag::Textured) {
        m_shader.bindTexture(*m_texture);
    }
    else{
        m_shader.setColor(m_color);
    }

    m_shader.draw(m_mesh);
}

