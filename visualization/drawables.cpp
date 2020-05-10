//
// Created by janos on 27.02.20.
//

#include "drawables.hpp"

#include <Magnum/Shaders/Flat.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/SceneGraph/Object.h>
#include <Magnum/SceneGraph/MatrixTransformation3D.h>
#include <Magnum/SceneGraph/Camera.h>


MeshDrawable::MeshDrawable(Object3D& obj, Mg::GL::Mesh& m, DrawableGroup* group, Mg::GL::Texture2D* t):
        Drawable(obj, group),
        mesh(m),
        texture(t)
{
}


FlatDrawable::FlatDrawable(Object3D& object, Mg::GL::Mesh& m, Mg::GL::AbstractShaderProgram& shader, DrawableGroup* group):
        MeshDrawable(object, m, group),
        shader(dynamic_cast<Magnum::Shaders::Flat3D&>(shader))
{
}

void FlatDrawable::draw(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera) {
    if (texture) shader.bindTexture(*texture);
    else shader.setColor(color);
    shader.setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix)
            .draw(mesh);
}

VertexColorDrawable::VertexColorDrawable(Object3D& object, Mg::GL::Mesh& m, Mg::GL::AbstractShaderProgram& shader, DrawableGroup* group):
        MeshDrawable(object, m, group),
        shader(&dynamic_cast<Mg::Shaders::VertexColor3D&>(shader))
{
}

void VertexColorDrawable::draw(const Magnum::Matrix4& transformationMatrix, Mg::SceneGraph::Camera3D& camera) {
    shader->setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix)
            .draw(mesh);
}



PhongDiffuseDrawable::PhongDiffuseDrawable(Object3D& object, Mg::GL::Mesh& m, Mg::GL::AbstractShaderProgram& shader, DrawableGroup* group):
        MeshDrawable(object, m, group),
        shader(dynamic_cast<Magnum::Shaders::Phong&>(shader))
{
}

void PhongDiffuseDrawable::draw(const Mg::Matrix4& transformationMatrix, Mg::SceneGraph::Camera3D& camera) {
    if(!texture) return;
    shader.bindDiffuseTexture(*texture)
            .setLightPosition({10.0f, 10.0f, 10.0f})
            .setTransformationMatrix(transformationMatrix)
            .setNormalMatrix(transformationMatrix.normalMatrix())
            .setProjectionMatrix(camera.projectionMatrix())
            .draw(mesh);
}



MeshVisualizerDrawable::MeshVisualizerDrawable(Object3D& object, Mg::GL::Mesh& m, Mg::GL::AbstractShaderProgram& s, DrawableGroup* group):
        MeshDrawable(object, m, group),
        shader(dynamic_cast<Mg::Shaders::MeshVisualizer3D&>(s))
{
}

void MeshVisualizerDrawable::draw(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera) {

    shader.setViewportSize(Mg::Vector2{Mg::GL::defaultFramebuffer.viewport().size()})
            .setTransformationMatrix(transformationMatrix)
            .setProjectionMatrix(camera.projectionMatrix())
            .setWireframeWidth(wireframeWidth)
            .setWireframeColor(wireframeColor)
            .setSmoothness(smoothness)
            .setColor(color)
            .draw(mesh);
}


FaceColorDrawable::FaceColorDrawable(Object3D& object, Mg::GL::Mesh& m, Mg::GL::AbstractShaderProgram& s, DrawableGroup* group):
        MeshDrawable(object, m, group),
        shader(dynamic_cast<Mg::Shaders::MeshVisualizer3D&>(s))
{
}

void FaceColorDrawable::draw(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera) {
    if(!texture) return;
    shader.setViewportSize(Mg::Vector2{Mg::GL::defaultFramebuffer.viewport().size()})
            .setTransformationMatrix(transformationMatrix)
            .setProjectionMatrix(camera.projectionMatrix())
            .setColorMapTransformation(offset, scale)
            .bindColorMapTexture(*texture)
            .draw(mesh);
}
