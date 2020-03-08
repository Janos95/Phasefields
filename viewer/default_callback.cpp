//
// Created by janos on 27.02.20.
//

#include "default_callback.hpp"
#include "object.hpp"

#include <Magnum/Shaders/Flat.h>

using namespace Magnum;
using namespace Corrade;

template<class... Ts> struct overloaded : Ts... { using Ts::operator()...; }; //TODO replace flat by variant

DefaultCallback::DefaultCallback(Object& object, Shaders::Flat3D& shader) :
        m_mesh(object.mesh),
        m_color(object.color),
        m_shader(shader)
{
    if(object.texture){
        m_texture = std::addressof(*object.texture);
    }
}

void DefaultCallback::operator()(const Matrix4& transformationMatrix, SceneGraph::Camera3D& camera){
    m_shader.setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix);

    if (m_texture) {
        m_shader.bindTexture(*m_texture);
    }
    else{
        m_shader.setColor(m_color);
    }

    m_mesh.draw(m_shader);
}
