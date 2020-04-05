//
// Created by janos on 27.02.20.
//

#pragma once

#include "drawable_data.hpp"

#include <Magnum/Math/Matrix4.h>
#include <Magnum/SceneGraph/Camera.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Color.h>
#include <Magnum/Shaders/Flat.h>
#include <Magnum/Shaders/VertexColor.h>
#include <Magnum/Shaders/Phong.h>
#include <Magnum/Shaders/MeshVisualizer.h>

class FlatCallback
{
public:

    explicit FlatCallback(TexturedDrawableData&, Magnum::Shaders::Flat3D&);

    void operator()(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera);

private:

    Magnum::GL::Mesh& m_mesh;
    Magnum::GL::Texture2D& m_texture;
    Magnum::Shaders::Flat3D& m_shader;
};

class VertexColorCallback
{
public:

    explicit VertexColorCallback(DrawableData& obj, Magnum::Shaders::VertexColor3D& shader): m_mesh(obj.mesh), m_shader(shader){}

    void operator()(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera){
        m_shader.setTransformationProjectionMatrix(camera.projectionMatrix() * transformationMatrix)
                .draw(m_mesh);
    }

private:
    Magnum::GL::Mesh& m_mesh;
    Magnum::Shaders::VertexColor3D& m_shader;
};

class PhongCallback
{
public:

    explicit PhongCallback(TexturedDrawableData& obj, Magnum::Shaders::Phong& shader):
        m_mesh(obj.mesh),
        m_shader(shader),
        m_texture(obj.texture)
        {
        }

    void operator()(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera){
        m_shader.bindDiffuseTexture(m_texture)
                .setLightPosition({5.0f, 5.0f, 7.0f})
                .setTransformationMatrix(transformationMatrix)
                .setNormalMatrix(transformationMatrix.normalMatrix())
                .setProjectionMatrix(camera.projectionMatrix())
                .draw(m_mesh);
    }

private:

    Magnum::GL::Mesh& m_mesh;
    Magnum::Shaders::Phong& m_shader;
    Magnum::GL::Texture2D& m_texture;
};

class MeshVisualizerCallback
{
public:

    explicit MeshVisualizerCallback(DrawableData& obj, Magnum::Shaders::MeshVisualizer3D& shader): m_mesh(obj.mesh), m_shader(shader){}

    void operator()(const Magnum::Matrix4& transformationMatrix, Magnum::SceneGraph::Camera3D& camera){
        m_shader.setColor(m_colorMesh)
                .setWireframeColor(m_colorWireframe)
                .setTransformationMatrix(transformationMatrix)
                .setProjectionMatrix(camera.projectionMatrix())
                .draw(m_mesh);
    }

private:

    Magnum::GL::Mesh& m_mesh;
    Magnum::Color4 m_colorMesh = Magnum::Color4::blue();
    Magnum::Color4 m_colorWireframe = Magnum::Color4::red();
    Magnum::Shaders::MeshVisualizer3D& m_shader;
};
