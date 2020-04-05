//
// Created by janos on 27.02.20.
//

#include "draw_callbacks.hpp"
#include "drawable_data.hpp"

#include <Magnum/Shaders/Flat.h>

using namespace Magnum;
using namespace Corrade;

FlatCallback::FlatCallback(TexturedDrawableData& object, Shaders::Flat3D& shader) :
        m_mesh(object.mesh),
        m_texture(object.texture),
        m_shader(shader)
{
}

void FlatCallback::operator()(const Matrix4& transformationMatrix, SceneGraph::Camera3D& camera){
    m_shader.setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix)
            .bindTexture(m_texture)
            .draw(m_mesh);
}

